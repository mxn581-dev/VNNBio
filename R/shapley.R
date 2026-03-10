# =============================================================================
# shapley.R -- Structure-Aware Shapley Pathway Attribution
#
# Axiomatically fair attribution of model predictions to biological pathways.
# Exploits VNN structure: pathway activations are independent given input,
# so coalition evaluation reduces to output-layer recomputation.
#
# Satisfies Shapley axioms: Efficiency, Symmetry, Linearity, Null Player.
# =============================================================================


# =============================================================================
# Activation functions (R-side mirrors of Julia definitions)
# =============================================================================

.activation_fn <- function(name) {
    switch(name,
        "tanh"    = tanh,
        "relu"    = function(x) pmax(x, 0),
        "sigmoid" = function(x) 1 / (1 + exp(-x)),
        "gelu"    = function(x) x * stats::pnorm(x),
        "swish"   = function(x) x / (1 + exp(-x)),
        stop("Unknown activation function: ", name, call. = FALSE)
    )
}

.sigmoid <- function(x) 1 / (1 + exp(-x))


# =============================================================================
# Parameter extraction from Julia
# =============================================================================

#' Extract All VNN Parameters to R
#'
#' Pulls weights, biases, and masks from the Julia backend into a
#' structured R list. After extraction, all downstream computation
#' (Shapley attribution, forward passes) is pure R.
#'
#' @param model A trained \code{VNNModel} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{layers}}{List of per-layer parameter lists, each with
#'       \code{weight}, \code{bias}, and optionally \code{mask}.}
#'     \item{\code{n_layers}}{Number of layers.}
#'     \item{\code{activation}}{Activation function name from the architecture.}
#'     \item{\code{task}}{Task type ("classification" or "regression").}
#'   }
#'
#' @keywords internal
.extractVNNParams <- function(model) {

    model_ref <- model@julia_model_ref
    if (is.null(model_ref)) {
        stop("Model has no Julia reference. Train the model first.",
             call. = FALSE)
    }

    if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
        stop("JuliaConnectoR required to extract parameters.",
             call. = FALSE)
    }

    raw <- JuliaConnectoR::juliaCall("_vnn_extract_vnn_params", model_ref)

    n_layers <- as.integer(raw[["n_layers"]])
    layers <- vector("list", n_layers)

    for (i in seq_len(n_layers)) {
        W <- as.matrix(raw[[paste0("layer_", i, "_weight")]])
        b <- as.numeric(raw[[paste0("layer_", i, "_bias")]])
        has_mask <- as.logical(raw[[paste0("layer_", i, "_has_mask")]])

        layer <- list(weight = W, bias = b, mask = NULL)
        if (has_mask) {
            layer$mask <- as.matrix(raw[[paste0("layer_", i, "_mask")]])
        }
        layers[[i]] <- layer
    }

    arch <- architecture(model)

    list(
        layers     = layers,
        n_layers   = n_layers,
        activation = activationFunction(arch),
        task       = taskType(model)
    )
}


# =============================================================================
# Pure R forward pass
# =============================================================================

#' R-side VNN Forward Pass
#'
#' Computes model prediction for a single sample using pre-extracted
#' parameters, with optional pathway masking for Shapley coalition
#' evaluation.
#'
#' @param x Numeric vector of length \code{n_genes} (one sample).
#' @param params Parameter list from \code{.extractVNNParams}.
#' @param active_pathways Integer vector of active pathway indices for
#'   the first masked layer. If \code{NULL}, all pathways are active.
#'
#' @return Numeric scalar prediction.
#' @keywords internal
.forwardPassR <- function(x, params, active_pathways = NULL) {
  
  act_fn <- .activation_fn(params$activation)
  h <- as.numeric(x)
  
  for (i in seq_along(params$layers)) {
    layer <- params$layers[[i]]
    W <- layer$weight
    b <- layer$bias
    
    ## Apply mask if present
    if (!is.null(layer$mask)) {
      W_masked <- W * layer$mask
    } else {
      W_masked <- W
    }
    
    ## z = W' * h + b
    if (nrow(W_masked) == length(h)) {
      z <- as.numeric(crossprod(W_masked, h)) + b
    } else {
      z <- as.numeric(W_masked %*% h) + b
    }
    
    ## Apply activation: hidden layers get the specified activation,
    ## output layer gets sigmoid (classification) or identity (regression)
    if (i < params$n_layers) {
      z <- act_fn(z)
      
      ## Zero out inactive pathways AFTER activation
      ## This ensures h_j = 0 for inactive pathways (not sigma(b_j))
      if (!is.null(active_pathways) && !is.null(layer$mask)) {
        inactive <- setdiff(seq_len(length(z)), active_pathways)
        if (length(inactive) > 0) {
          z[inactive] <- 0
        }
      }
    } else {
      if (params$task == "classification") {
        z <- .sigmoid(z)
      }
    }
    
    h <- z
  }
  
  h
}


# =============================================================================
# Structure-aware precomputation
# =============================================================================

#' Precompute Per-Pathway Activations
#'
#' For a single-layer VNN, each pathway's hidden activation depends only
#' on its own genes (by construction of the mask). This function computes
#' all pathway activations independently, enabling O(1) coalition
#' evaluation at the output layer.
#'
#' @param x Numeric vector of length n_genes.
#' @param params Parameter list from .extractVNNParams.
#'
#' @return List with:
#'   \describe{
#'     \item{\code{h}}{Numeric vector of pathway activations.}
#'     \item{\code{W_out}}{Output layer weight matrix.}
#'     \item{\code{b_out}}{Output layer bias vector.}
#'     \item{\code{task}}{Task type.}
#'   }
#' @keywords internal
.precomputePathwayActivations <- function(x, params) {

    act_fn <- .activation_fn(params$activation)

    ## Layer 1: masked hidden layer
    l1 <- params$layers[[1]]
    W_masked <- l1$weight * l1$mask
    h <- act_fn(as.numeric(crossprod(W_masked, as.numeric(x))) + l1$bias)

    ## Output layer (last layer, no mask)
    l_out <- params$layers[[params$n_layers]]

    list(
        h     = h,
        W_out = l_out$weight,
        b_out = l_out$bias,
        task  = params$task
    )
}


#' Evaluate Coalition via Precomputed Activations
#'
#' Given precomputed pathway activations, compute the model output when
#' only pathways in the coalition are active. Inactive pathways contribute
#' zero to the output sum.
#'
#' @param precomp Output of \code{.precomputePathwayActivations}.
#' @param coalition Integer vector of active pathway indices.
#'
#' @return Numeric scalar prediction.
#' @keywords internal
.evalCoalition <- function(precomp, coalition) {
    ## h_active: only include coalition members
    h_masked <- rep(0, length(precomp$h))
    h_masked[coalition] <- precomp$h[coalition]

    ## Output layer: y = f(W_out' * h_masked + b_out)
    w_vec <- as.numeric(precomp$W_out)
    logit <- sum(w_vec * h_masked) + as.numeric(precomp$b_out)

    if (precomp$task == "classification") {
        .sigmoid(logit)
    } else {
        logit
    }
}


# =============================================================================
# Main Shapley computation
# =============================================================================

#' Shapley Pathway Attribution for Visible Neural Networks
#'
#' Computes axiomatically fair Shapley values attributing model predictions
#' to biological pathways. Exploits VNN structure for efficient computation:
#' since pathway activations are independent given the input, coalition
#' evaluation reduces to output-layer recomputation (O(1) per evaluation
#' instead of O(n_genes * n_pathways)).
#'
#' @details
#' ## Shapley Value Definition
#'
#' For pathway \eqn{j} in a model with pathway set \eqn{P}:
#'
#' \deqn{\phi_j = \frac{1}{|P|!} \sum_{\pi} \left[ v(S_j^{\pi} \cup \{j\})
#'   - v(S_j^{\pi}) \right]}
#'
#' where \eqn{\pi} ranges over all orderings of \eqn{P} and
#' \eqn{S_j^{\pi}} is the set of pathways preceding \eqn{j} in ordering
#' \eqn{\pi}.
#'
#' ## Structure-Aware Speedup
#'
#' For a single-layer VNN (genes → pathways → output), each pathway's
#' activation \eqn{h_j = \sigma(\sum_i W_{ij} M_{ij} x_i + b_j)} depends
#' only on genes in that pathway. We precompute all \eqn{h_j} once, then
#' evaluate any coalition \eqn{v(S)} as
#' \eqn{f(\sum_{j \in S} w_j^{out} h_j + b^{out})} in O(1).
#'
#' Total cost: O(M * P * N) where M = permutations, P = pathways,
#' N = samples. For typical values (M=200, P=50, N=663), this completes
#' in seconds.
#'
#' ## Shapley Axioms Satisfied
#'
#' \describe{
#'   \item{Efficiency}{\eqn{\sum_j \phi_j = f(x) - f(\emptyset)}.
#'     Attributions sum to the difference between the full prediction
#'     and the empty-coalition baseline.}
#'   \item{Symmetry}{If two pathways contribute identically to every
#'     coalition, they receive equal attribution.}
#'   \item{Null Player}{A pathway whose activation is zero (or whose
#'     output weight is zero) receives zero attribution.}
#'   \item{Linearity}{Shapley values of a sum of games equal the sum
#'     of Shapley values.}
#' }
#'
#' @param model A trained \code{VNNModel} object.
#' @param newdata A \code{SummarizedExperiment} or numeric matrix
#'   (\code{[samples x genes]}) for which to compute attributions.
#' @param mode One of \code{"global"} (average across samples),
#'   \code{"local"} (per-sample attributions), or \code{"both"}.
#' @param n_perm Integer. Number of permutation samples for Shapley
#'   approximation. Default 200. More permutations = more precise but
#'   slower.
#' @param assay_name Character. Assay to use from SE. Default first.
#' @param params Optional pre-extracted parameter list from
#'   \code{.extractVNNParams}. If provided, Julia is not called.
#'   Enables fully Julia-free operation with saved parameters.
#' @param seed Integer. Random seed. Default 42L.
#' @param verbose Logical. Print progress? Default TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{global}}{Named numeric vector of mean absolute Shapley
#'       values per pathway (if \code{mode} is "global" or "both").}
#'     \item{\code{local}}{Matrix of Shapley values
#'       \code{[n_samples x n_pathways]} (if \code{mode} is "local" or
#'       "both"). Positive = pushes prediction toward 1, negative = toward 0.}
#'     \item{\code{predictions}}{Numeric vector of full-model predictions.}
#'     \item{\code{baseline}}{Scalar prediction with no pathways active.}
#'     \item{\code{pathway_names}}{Character vector of pathway names.}
#'     \item{\code{n_perm}}{Number of permutations used.}
#'   }
#'
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay assayNames
#'
#' @export
#' @examples
#' \donttest{
#' # After training:
#' # shap <- shapleyPathwayAttribution(model, se, mode = "both")
#' # shap$global          # mean |SHAP| per pathway
#' # shap$local[1, ]      # patient 1's pathway attributions
#' # plotShapleyGlobal(shap)
#' # plotShapleyLocal(shap, sample_idx = 1)
#' }
shapleyPathwayAttribution <- function(model,
                                       newdata,
                                       mode = c("global", "local", "both"),
                                       n_perm = 200L,
                                       assay_name = NULL,
                                       params = NULL,
                                       seed = 42L,
                                       verbose = TRUE) {

    mode <- match.arg(mode)
    n_perm <- as.integer(n_perm)
    set.seed(seed)

    ## ---- input validation ---------------------------------------------------
    stopifnot(is(model, "VNNModel"))

    if (n_perm < 10L) {
        stop("n_perm must be >= 10 for meaningful Shapley estimates.",
             call. = FALSE)
    }

    ## ---- extract parameters -------------------------------------------------
    if (is.null(params)) {
        params <- .extractVNNParams(model)
    }

    ## Verify we have at least 2 layers (1 masked + 1 output)
    if (params$n_layers < 2) {
        stop("Shapley attribution requires at least 2 layers ",
             "(masked hidden + output). Found ", params$n_layers, ".",
             call. = FALSE)
    }

    ## Check structure-aware eligibility: first layer must have mask
    use_fast <- !is.null(params$layers[[1]]$mask)
    if (!use_fast) {
        stop("First layer has no mask. Structure-aware Shapley requires ",
             "at least one masked layer.", call. = FALSE)
    }

    ## ---- extract expression matrix ------------------------------------------
    if (is(newdata, "SummarizedExperiment")) {
        if (is.null(assay_name)) {
            assay_name <- SummarizedExperiment::assayNames(newdata)[1]
        }
        X <- t(as.matrix(SummarizedExperiment::assay(newdata, assay_name)))
    } else if (is.matrix(newdata)) {
        X <- newdata
    } else {
        stop("'newdata' must be a SummarizedExperiment or [samples x genes] ",
             "matrix.", call. = FALSE)
    }

    n_samples <- nrow(X)
    n_pathways <- ncol(params$layers[[1]]$mask)

    ## Get pathway names from architecture mask
    arch <- architecture(model)
    masks <- layerMasks(arch)
    pw_names <- colnames(masks[[1]])
    if (is.null(pw_names)) {
        pw_names <- paste0("pathway_", seq_len(n_pathways))
    }

    ## ---- compute Shapley values ---------------------------------------------
    shap_matrix <- matrix(0, nrow = n_samples, ncol = n_pathways)
    colnames(shap_matrix) <- pw_names
    predictions <- numeric(n_samples)
    all_pathways <- seq_len(n_pathways)

    if (verbose) cat("Computing Shapley values for", n_samples, "samples...\n")

    for (s in seq_len(n_samples)) {
        if (verbose && s %% 10 == 0) {
            cat(sprintf("  Sample %d/%d\r", s, n_samples))
        }

        x_s <- as.numeric(X[s, ])

        ## Precompute pathway activations (structure-aware speedup)
        precomp <- .precomputePathwayActivations(x_s, params)

        ## Full prediction (all pathways active)
        predictions[s] <- .evalCoalition(precomp, all_pathways)

        ## Permutation-based Shapley approximation
        marginals <- matrix(0, nrow = n_perm, ncol = n_pathways)

        for (m in seq_len(n_perm)) {
            ## Random ordering of pathways
            perm <- sample(n_pathways)

            ## Walk through the ordering, tracking the running coalition
            ## and computing marginal contributions
            running_val <- .evalCoalition(precomp, integer(0))

            for (pos in seq_along(perm)) {
                j <- perm[pos]
                coalition <- perm[seq_len(pos)]
                new_val <- .evalCoalition(precomp, coalition)
                marginals[m, j] <- new_val - running_val
                running_val <- new_val
            }
        }

        ## Average marginal contributions across permutations
        shap_matrix[s, ] <- colMeans(marginals)
    }

    if (verbose) cat("\n")

    ## ---- baseline (empty coalition) -----------------------------------------
    ## Use first sample's precomp structure (output layer params are shared)
    precomp_ref <- .precomputePathwayActivations(as.numeric(X[1, ]), params)
    baseline <- .evalCoalition(precomp_ref, integer(0))

    ## ---- assemble results ---------------------------------------------------
    result <- list(
        pathway_names = pw_names,
        predictions   = predictions,
        baseline      = baseline,
        n_perm        = n_perm
    )

    if (mode %in% c("local", "both")) {
        result$local <- shap_matrix
    }

    if (mode %in% c("global", "both")) {
        ## Global importance = mean |SHAP| across samples
        global <- colMeans(abs(shap_matrix))
        result$global <- sort(global, decreasing = TRUE)
    }

    if (verbose) {
        cat(sprintf("Shapley attribution complete.\n"))
        cat(sprintf("  Samples: %d | Pathways: %d | Permutations: %d\n",
                    n_samples, n_pathways, n_perm))

        ## Verify efficiency axiom (attributions sum to prediction - baseline)
        efficiency_check <- mean(abs(
            rowSums(shap_matrix) - (predictions - baseline)
        ))
        cat(sprintf("  Efficiency check (mean |residual|): %.6f\n",
                    efficiency_check))
        if (efficiency_check > 0.01) {
            warning("Efficiency residual > 0.01. Consider increasing n_perm.",
                    call. = FALSE)
        }
    }

    result
}


# =============================================================================
# Visualization: Global Shapley importance
# =============================================================================

#' Plot Global Shapley Pathway Importance
#'
#' Bar chart of mean |SHAP| per pathway, showing which pathways have the
#' largest average impact on predictions across all samples.
#'
#' @param shap_result Result list from \code{shapleyPathwayAttribution}
#'   (must contain \code{global}).
#' @param top_n Integer. Number of top pathways to display. Default 20.
#'
#' @return A \code{ggplot} object.
#' @export
plotShapleyGlobal <- function(shap_result, top_n = 20L) {

    if (is.null(shap_result$global)) {
        stop("No global Shapley values found. ",
             "Use mode = 'global' or 'both' in shapleyPathwayAttribution().",
             call. = FALSE)
    }

    global <- shap_result$global
    n_show <- min(top_n, length(global))
    top <- global[seq_len(n_show)]

    plot_df <- data.frame(
        pathway    = names(top),
        importance = unname(top),
        stringsAsFactors = FALSE
    )
    plot_df$pathway <- factor(plot_df$pathway,
                              levels = rev(plot_df$pathway))

    ggplot(plot_df, aes(x = .data$pathway, y = .data$importance,
                        fill = .data$importance)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7",
                             high = "#b2182b",
                             midpoint = max(plot_df$importance) / 2) +
        labs(
            title = "Shapley Pathway Attribution (Global)",
            subtitle = paste0("Mean |SHAP| across ", nrow(shap_result$local),
                              " samples"),
            x = NULL,
            y = "Mean |Shapley Value|"
        ) +
        theme_minimal(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
}


# =============================================================================
# Visualization: Local (per-sample) Shapley
# =============================================================================

#' Plot Local Shapley Attribution for a Single Sample
#'
#' Horizontal bar chart showing each pathway's signed Shapley contribution
#' to a single sample's prediction. Positive values push the prediction
#' toward class 1, negative toward class 0.
#'
#' @param shap_result Result list from \code{shapleyPathwayAttribution}
#'   (must contain \code{local}).
#' @param sample_idx Integer. Which sample to display (row index).
#' @param top_n Integer. Show top N pathways by |SHAP|. Default 15.
#' @param sample_label Optional character label for the sample.
#'
#' @return A \code{ggplot} object.
#' @export
plotShapleyLocal <- function(shap_result,
                              sample_idx = 1L,
                              top_n = 15L,
                              sample_label = NULL) {

    if (is.null(shap_result$local)) {
        stop("No local Shapley values found. ",
             "Use mode = 'local' or 'both' in shapleyPathwayAttribution().",
             call. = FALSE)
    }

    shap_row <- shap_result$local[sample_idx, ]
    pred <- shap_result$predictions[sample_idx]

    ## Select top pathways by absolute value
    ord <- order(abs(shap_row), decreasing = TRUE)
    n_show <- min(top_n, length(shap_row))
    top_idx <- ord[seq_len(n_show)]

    plot_df <- data.frame(
        pathway = names(shap_row)[top_idx],
        shap    = unname(shap_row[top_idx]),
        stringsAsFactors = FALSE
    )
    plot_df$direction <- ifelse(plot_df$shap >= 0, "Positive", "Negative")

    ## Order by signed value for visual clarity
    plot_df <- plot_df[order(plot_df$shap), ]
    plot_df$pathway <- factor(plot_df$pathway, levels = plot_df$pathway)

    label <- if (!is.null(sample_label)) {
        sample_label
    } else {
        paste0("Sample ", sample_idx)
    }

    ggplot(plot_df, aes(x = .data$pathway, y = .data$shap,
                        fill = .data$direction)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        geom_hline(yintercept = 0, color = "gray40", linewidth = 0.3) +
        ggplot2::scale_fill_manual(
            values = c("Positive" = "#b2182b", "Negative" = "#2166ac")
        ) +
        labs(
            title = paste0("Shapley Attribution: ", label),
            subtitle = sprintf("Prediction: %.3f | Baseline: %.3f",
                               pred, shap_result$baseline),
            x = NULL,
            y = "Shapley Value (contribution to prediction)"
        ) +
        theme_minimal(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
}


# =============================================================================
# Visualization: Beeswarm / summary plot (all samples, all pathways)
# =============================================================================

#' Plot Shapley Summary (Beeswarm-style)
#'
#' Shows the distribution of Shapley values across all samples for the top
#' pathways, combining global importance ranking with local value spread.
#' Similar to SHAP summary plots in Python.
#'
#' @param shap_result Result list from \code{shapleyPathwayAttribution}
#'   (must contain \code{local}).
#' @param top_n Integer. Number of top pathways. Default 15.
#'
#' @return A \code{ggplot} object.
#' @export
plotShapleySummary <- function(shap_result, top_n = 15L) {

    if (is.null(shap_result$local)) {
        stop("No local Shapley values found. ",
             "Use mode = 'local' or 'both' in shapleyPathwayAttribution().",
             call. = FALSE)
    }

    ## Rank by mean |SHAP|
    mean_abs <- colMeans(abs(shap_result$local))
    ord <- order(mean_abs, decreasing = TRUE)
    n_show <- min(top_n, length(mean_abs))
    top_idx <- ord[seq_len(n_show)]

    ## Build long-form data.frame
    rows <- lapply(top_idx, function(j) {
        data.frame(
            pathway = shap_result$pathway_names[j],
            shap    = shap_result$local[, j],
            stringsAsFactors = FALSE
        )
    })
    plot_df <- do.call(rbind, rows)

    ## Order pathways by mean |SHAP| (bottom = most important for coord_flip)
    pw_order <- shap_result$pathway_names[top_idx]
    plot_df$pathway <- factor(plot_df$pathway, levels = rev(pw_order))

    ggplot(plot_df, aes(x = .data$pathway, y = .data$shap)) +
        ggplot2::geom_jitter(
            width = 0.2, height = 0, alpha = 0.3, size = 0.8,
            color = "#b2182b"
        ) +
        ggplot2::stat_summary(
            fun = median, geom = "point",
            shape = 18, size = 2.5, color = "black"
        ) +
        coord_flip() +
        geom_hline(yintercept = 0, color = "gray40", linewidth = 0.3) +
        labs(
            title = "Shapley Value Distribution by Pathway",
            subtitle = paste0("Top ", n_show, " pathways | ",
                              "Diamonds = median"),
            x = NULL,
            y = "Shapley Value"
        ) +
        theme_minimal(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
}
