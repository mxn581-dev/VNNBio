# =============================================================================
# pretrained.R -- ARCHS4 Pretrained Feature Extraction (Mode 1: Frozen)
#
# Pure R, no Julia dependency. Uses ARCHS4-pretrained weights as a fixed
# feature extractor: genes → pathway activations via masked hidden layer,
# then glmnet for classification on the pathway-level representation.
#
# Mode 1 pipeline:
#   1. loadPretrained("archs4_hallmark_v1")  →  PretrainedVNN object
#   2. classifyFrozen(pretrained, se, "label")  →  predictions + AUC
# =============================================================================


# ── S4 Class ─────────────────────────────────────────────────────

#' PretrainedVNN: Frozen ARCHS4-Pretrained Feature Extractor
#'
#' Stores pretrained hidden-layer weights, bias, and mask from ARCHS4
#' multi-task training. Used as a fixed feature extractor that maps
#' gene expression to pathway-level activations without any Julia
#' dependency.
#'
#' @slot W_hidden Dense numeric matrix \code{[n_genes x n_pathways]}.
#'   Pretrained weights from ARCHS4 multi-task learning.
#' @slot b_hidden Numeric vector of length \code{n_pathways}. Pretrained
#'   bias terms.
#' @slot mask Numeric matrix \code{[n_genes x n_pathways]}. Binary mask
#'   constraining gene-pathway connectivity (from MSigDB Hallmark).
#' @slot gene_info A \code{data.frame} with columns \code{archs4_index},
#'   \code{symbol}, \code{ensembl}. Row order matches \code{W_hidden}.
#' @slot pathway_names Character vector of pathway names.
#' @slot metadata Named list of provenance info (source, accuracy, etc.).
#'
#' @section Constructors:
#' Use \code{\link{loadPretrained}} rather than calling
#' \code{new("PretrainedVNN", ...)} directly.
#'
#' @name PretrainedVNN-class
#' @rdname PretrainedVNN-class
#' @exportClass PretrainedVNN
setClass("PretrainedVNN",
    slots = list(
        W_hidden      = "matrix",
        b_hidden      = "numeric",
        mask          = "matrix",
        gene_info     = "data.frame",
        pathway_names = "character",
        metadata      = "list"
    )
)

setValidity("PretrainedVNN", function(object) {
    errors <- character()
    if (nrow(object@W_hidden) != nrow(object@mask))
        errors <- c(errors, "W_hidden and mask must have same number of rows")
    if (ncol(object@W_hidden) != ncol(object@mask))
        errors <- c(errors, "W_hidden and mask must have same number of columns")
    if (ncol(object@W_hidden) != length(object@b_hidden))
        errors <- c(errors, "ncol(W_hidden) must equal length(b_hidden)")
    if (nrow(object@W_hidden) != nrow(object@gene_info))
        errors <- c(errors, "nrow(W_hidden) must equal nrow(gene_info)")
    if (ncol(object@W_hidden) != length(object@pathway_names))
        errors <- c(errors, "ncol(W_hidden) must equal length(pathway_names)")

    req_cols <- c("archs4_index", "symbol", "ensembl")
    if (!all(req_cols %in% colnames(object@gene_info)))
        errors <- c(errors, paste0("gene_info must have columns: ",
                                   paste(req_cols, collapse = ", ")))

    if (length(errors) == 0) TRUE else errors
})

#' @rdname PretrainedVNN-class
#' @aliases show,PretrainedVNN-method
setMethod("show", "PretrainedVNN", function(object) {
    cat("PretrainedVNN object\n")
    cat("  Weights ID: ", object@metadata$weights_id %||% "unknown", "\n")
    cat("  Source:     ", object@metadata$source %||% "unknown", "\n")
    cat("  Genes:      ", nrow(object@W_hidden), "\n")
    cat("  Pathways:   ", ncol(object@W_hidden), "\n")
    cat("  Activation:  tanh (frozen)\n")
    if (!is.null(object@metadata$tissue_accuracy))
        cat("  Pretrain tissue acc:  ",
            round(object@metadata$tissue_accuracy * 100, 1), "%\n")
    if (!is.null(object@metadata$disease_accuracy))
        cat("  Pretrain disease acc: ",
            round(object@metadata$disease_accuracy * 100, 1), "%\n")
    n_show <- min(5, length(object@pathway_names))
    cat("  Pathways: ",
        paste(head(object@pathway_names, n_show), collapse = ", "),
        if (length(object@pathway_names) > n_show) ", ..." else "", "\n")
})


# ── Load pretrained weights ──────────────────────────────────────

#' Load Pretrained ARCHS4 Weights
#'
#' Loads a pretrained weight bundle into a \code{PretrainedVNN} object.
#' By default, looks for bundled weights in the package's
#' \code{inst/extdata/} directory. Alternatively, a local file path
#' can be provided for weights stored outside the package (e.g.,
#' distributed via GitHub release).
#'
#' @param weights_id Character. Identifier for the weight bundle.
#'   Default \code{"archs4_hallmark_v1"}.
#' @param path Optional file path to a \code{.rds} weight bundle.
#'   If \code{NULL}, the function searches \code{inst/extdata/} for
#'   \code{<weights_id>.rds}.
#'
#' @return A \code{PretrainedVNN} object.
#'
#' @details
#' The pretrained weights were learned by training a VNN on the full
#' ARCHS4 compendium (~500K human RNA-seq samples) with tissue type
#' and disease state as multi-task labels, using MSigDB Hallmark
#' pathways as the architecture constraint.
#'
#' The returned object acts as a frozen feature extractor: pass it to
#' \code{\link{classifyFrozen}} to classify new datasets using the
#' pretrained pathway representation, or to
#' \code{\link{trainVNN}} via the \code{pretrained} argument for
#' fine-tuned transfer learning (Mode 2).
#'
#' @export
#' @examples
#' \donttest{
#' pt <- loadPretrained("archs4_hallmark_v1")
#' pt
#' }
loadPretrained <- function(weights_id = "archs4_hallmark_v1",
                           path = NULL) {

    if (is.null(path)) {
        filename <- paste0(weights_id, ".rds")
        path <- system.file("extdata", filename,
                            package = "VNNBio", mustWork = FALSE)
        if (!nzchar(path)) {
            stop(sprintf(
                paste0("Pretrained weights '%s' not found in package.\n",
                       "  Expected: inst/extdata/%s\n",
                       "  Provide a local path via the 'path' argument,\n",
                       "  or run data-raw/build_pretrained_weights.R first."),
                weights_id, filename
            ), call. = FALSE)
        }
    }

    if (!file.exists(path)) {
        stop("File not found: ", path, call. = FALSE)
    }

    bundle <- readRDS(path)

    ## Validate structure
    required_keys <- c("W_hidden", "b_hidden", "mask", "gene_info",
                       "pathway_names", "metadata")
    missing_keys <- setdiff(required_keys, names(bundle))
    if (length(missing_keys) > 0) {
        stop("Pretrained bundle is missing keys: ",
             paste(missing_keys, collapse = ", "),
             call. = FALSE)
    }

    obj <- new("PretrainedVNN",
        W_hidden      = bundle$W_hidden,
        b_hidden      = bundle$b_hidden,
        mask          = bundle$mask,
        gene_info     = bundle$gene_info,
        pathway_names = bundle$pathway_names,
        metadata      = bundle$metadata
    )

    validObject(obj)
    message(sprintf("Loaded pretrained weights: %s (%d genes, %d pathways)",
                    weights_id, nrow(obj@W_hidden), ncol(obj@W_hidden)))
    obj
}


# ── Gene alignment helper ────────────────────────────────────────

#' Align User Gene Space to Pretrained Gene Space
#'
#' Identifies the intersection between the user's expression features
#' and the pretrained model's gene list, then subsets both the
#' expression matrix and the pretrained weight/mask matrices to this
#' shared space. Handles both gene symbol and Ensembl ID matching.
#'
#' @param se A \code{SummarizedExperiment}.
#' @param pretrained A \code{PretrainedVNN} object.
#' @param assay_name Character. Which assay to use.
#' @param gene_id_type One of \code{"symbol"}, \code{"ensembl"}, or
#'   \code{"auto"} (default). Auto-detection checks if SE rownames
#'   start with "ENSG".
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{X}}{Numeric matrix \code{[samples x aligned_genes]}.}
#'     \item{\code{W}}{Pretrained weights \code{[aligned_genes x pathways]}.}
#'     \item{\code{b}}{Pretrained bias \code{[pathways]}.}
#'     \item{\code{mask}}{Binary mask \code{[aligned_genes x pathways]}.}
#'     \item{\code{n_matched}}{Number of genes successfully aligned.}
#'     \item{\code{n_user}}{Total genes in user SE.}
#'     \item{\code{gene_id_type}}{Detected or specified ID type.}
#'   }
#'
#' @keywords internal
.alignGenesToPretrained <- function(se, pretrained, assay_name,
                                    gene_id_type = "auto") {

    user_genes <- rownames(se)
    if (is.null(user_genes) || length(user_genes) == 0) {
        stop("SummarizedExperiment must have rownames (gene IDs).",
             call. = FALSE)
    }

    ## Auto-detect ID type
    if (gene_id_type == "auto") {
        gene_id_type <- if (any(grepl("^ENSG", head(user_genes, 20)))) {
            "ensembl"
        } else {
            "symbol"
        }
    }

    ## Select the matching column from gene_info
    pt_col <- switch(gene_id_type,
        symbol  = "symbol",
        ensembl = "ensembl",
        stop("gene_id_type must be 'symbol', 'ensembl', or 'auto'.",
             call. = FALSE)
    )

    pt_ids <- pretrained@gene_info[[pt_col]]

    ## Strip Ensembl version suffixes from user genes if needed
    user_genes_clean <- user_genes
    if (gene_id_type == "ensembl") {
        user_genes_clean <- sub("\\.\\d+$", "", user_genes_clean)
    }

    ## Find intersection (preserving pretrained order for weight alignment)
    shared <- intersect(pt_ids, user_genes_clean)

    if (length(shared) == 0) {
        stop(sprintf(
            paste0("No gene overlap between SE (%s IDs) and pretrained ",
                   "model (%s column).\n",
                   "  First 5 SE genes: %s\n",
                   "  First 5 pretrained: %s\n",
                   "  Check gene_id_type or strip Ensembl version suffixes."),
            gene_id_type, pt_col,
            paste(head(user_genes_clean, 5), collapse = ", "),
            paste(head(pt_ids, 5), collapse = ", ")
        ), call. = FALSE)
    }

    ## Indices into pretrained space
    pt_idx <- match(shared, pt_ids)

    ## Indices into user SE (map back to original rownames)
    user_lookup <- setNames(seq_along(user_genes_clean), user_genes_clean)
    user_idx <- unname(user_lookup[shared])

    ## Subset expression matrix: SE is [genes x samples], transpose
    X_full <- as.matrix(SummarizedExperiment::assay(se, assay_name))
    X <- t(X_full[user_idx, , drop = FALSE])  # [samples x aligned_genes]

    ## Subset pretrained weights and mask
    W_aligned    <- pretrained@W_hidden[pt_idx, , drop = FALSE]
    mask_aligned <- pretrained@mask[pt_idx, , drop = FALSE]

    ## Report coverage
    coverage <- length(shared) / length(pt_ids) * 100
    message(sprintf(
        "Gene alignment: %d / %d pretrained genes matched (%.1f%% coverage)",
        length(shared), length(pt_ids), coverage
    ))

    if (coverage < 50) {
        warning(sprintf(
            paste0("Low gene coverage (%.1f%%). Pathway activations may be ",
                   "unreliable. Consider using %s IDs."),
            coverage,
            if (gene_id_type == "symbol") "ensembl" else "symbol"
        ), call. = FALSE)
    }

    list(
        X            = X,
        W            = W_aligned,
        b            = pretrained@b_hidden,
        mask         = mask_aligned,
        n_matched    = length(shared),
        n_user       = length(user_genes),
        gene_id_type = gene_id_type
    )
}


# ── Mode 1: Frozen classification ────────────────────────────────

#' Classify Using Frozen Pretrained Pathway Features
#'
#' Uses ARCHS4-pretrained weights as a fixed feature extractor to
#' transform gene expression into pathway-level activations, then
#' fits a penalized logistic regression (glmnet) for classification.
#' Entirely pure R — no Julia required.
#'
#' @details
#' ## Pipeline
#'
#' 1. **Gene alignment**: Intersect user genes with pretrained gene list.
#' 2. **Forward pass**: \eqn{H = \tanh(X \cdot (W \odot M) + b)} where
#'    \eqn{X} is \code{[samples x genes]}, \eqn{W \odot M} is the
#'    masked pretrained weight matrix, and \eqn{b} is the pretrained
#'    bias. This produces \eqn{H} of shape \code{[samples x pathways]}.
#' 3. **Classification**: \code{cv.glmnet} with \code{family="binomial"}
#'    on \eqn{H} predicts the binary label. Uses leave-one-out
#'    cross-validation by default for small datasets.
#' 4. **Pathway importance**: glmnet coefficients directly quantify each
#'    pathway's contribution to classification.
#'
#' ## Why This Works
#'
#' The ARCHS4 pretraining learned gene-to-pathway weight patterns across
#' ~500K samples, capturing general transcriptomic structure. The frozen
#' hidden layer acts as a biologically-informed dimensionality reduction
#' (e.g., 20,000 genes → 50 Hallmark pathways). The downstream glmnet
#' then only needs to learn a simple linear boundary in this compact,
#' pre-structured space — making it effective even with very small
#' sample sizes (n < 100).
#'
#' @param pretrained A \code{PretrainedVNN} object from
#'   \code{\link{loadPretrained}}.
#' @param se A \code{SummarizedExperiment} with gene expression.
#' @param label_col Character. Column in \code{colData(se)} containing
#'   the binary class label.
#' @param assay_name Character. Which assay to use. Default: first
#'   available.
#' @param gene_id_type One of \code{"symbol"}, \code{"ensembl"}, or
#'   \code{"auto"} (default). Controls how SE rownames are matched to
#'   pretrained genes.
#' @param nfolds Integer. Number of CV folds for glmnet. Default
#'   \code{NULL} = leave-one-out (LOOCV). Set to e.g. 10 for larger
#'   datasets.
#' @param alpha Numeric. Elastic net mixing parameter for glmnet
#'   (1 = lasso, 0 = ridge). Default 0.5.
#' @param seed Integer. Random seed. Default 42L.
#' @param verbose Logical. Print progress? Default TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{predictions}}{Numeric vector of predicted
#'       probabilities for each sample.}
#'     \item{\code{labels}}{Factor of true labels.}
#'     \item{\code{auc}}{AUC from cross-validated predictions.}
#'     \item{\code{pathway_activations}}{Matrix \code{[samples x pathways]}
#'       of hidden-layer activations \eqn{H}.}
#'     \item{\code{pathway_importance}}{Named numeric vector of absolute
#'       glmnet coefficients, sorted descending.}
#'     \item{\code{glmnet_fit}}{The \code{cv.glmnet} object for
#'       further inspection.}
#'     \item{\code{gene_coverage}}{Fraction of pretrained genes matched.}
#'     \item{\code{lambda}}{Selected lambda value (\code{lambda.min}).}
#'   }
#'
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom stats predict
#'
#' @export
#' @examples
#' \donttest{
#' pt <- loadPretrained("archs4_hallmark_v1")
#' result <- classifyFrozen(pt, se, label_col = "condition")
#' result$auc
#' head(result$pathway_importance)
#' }
classifyFrozen <- function(pretrained,
                           se,
                           label_col,
                           assay_name = NULL,
                           gene_id_type = "auto",
                           nfolds = NULL,
                           alpha = 0.5,
                           seed = 42L,
                           verbose = TRUE) {

    ## ---- dependency check ----------------------------------------------------
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required for classifyFrozen().\n",
             "Install with: install.packages('glmnet')",
             call. = FALSE)
    }

    ## ---- input validation ----------------------------------------------------
    stopifnot(is(pretrained, "PretrainedVNN"))
    stopifnot(is(se, "SummarizedExperiment"))

    if (!label_col %in% colnames(SummarizedExperiment::colData(se))) {
        stop("label_col '", label_col, "' not found in colData(se).",
             call. = FALSE)
    }

    if (is.null(assay_name)) {
        assay_name <- SummarizedExperiment::assayNames(se)[1]
    }

    ## ---- extract and encode labels -------------------------------------------
    y_raw <- SummarizedExperiment::colData(se)[[label_col]]
    y_factor <- as.factor(y_raw)
    levels_y <- levels(y_factor)

    if (length(levels_y) != 2) {
        stop(sprintf(
            "classifyFrozen requires exactly 2 classes. Found %d: %s",
            length(levels_y), paste(levels_y, collapse = ", ")
        ), call. = FALSE)
    }

    y_binary <- as.integer(y_factor) - 1L  # 0/1

    if (verbose) {
        cat(sprintf("Labels: '%s' (n=%d) vs '%s' (n=%d)\n",
                    levels_y[1], sum(y_binary == 0),
                    levels_y[2], sum(y_binary == 1)))
    }

    ## ---- align genes ---------------------------------------------------------
    aligned <- .alignGenesToPretrained(se, pretrained, assay_name,
                                       gene_id_type)
    X <- aligned$X  # [samples x aligned_genes]

    ## ---- frozen forward pass -------------------------------------------------
    ## H = tanh(X %*% (W * mask) + b)
    W_masked <- aligned$W * aligned$mask  # [aligned_genes x pathways]

    # Broadcast bias: each row of (X %*% W_masked) gets b added
    H_raw <- X %*% W_masked  # [samples x pathways]
    H <- tanh(sweep(H_raw, 2, aligned$b, "+"))

    colnames(H) <- pretrained@pathway_names

    if (verbose) {
        cat(sprintf("Pathway activations: [%d samples x %d pathways]\n",
                    nrow(H), ncol(H)))
    }

    ## ---- glmnet classification -----------------------------------------------
    set.seed(seed)

    # LOOCV for small datasets, k-fold for larger
    if (is.null(nfolds)) {
        nfolds <- nrow(H)  # LOOCV
    }
    nfolds <- min(nfolds, nrow(H))

    if (verbose) {
        fold_label <- if (nfolds == nrow(H)) "LOOCV" else paste0(nfolds, "-fold CV")
        cat(sprintf("Fitting glmnet (alpha=%.2f, %s)...\n", alpha, fold_label))
    }

    cv_fit <- glmnet::cv.glmnet(
        x      = H,
        y      = y_binary,
        family  = "binomial",
        alpha   = alpha,
        nfolds  = nfolds,
        type.measure = "auc",
        standardize  = TRUE
    )

    ## ---- extract predictions -------------------------------------------------
    preds <- as.numeric(
        stats::predict(cv_fit, newx = H, s = "lambda.min", type = "response")
    )

    ## ---- compute AUC ---------------------------------------------------------
    # Use the cv.glmnet's own cross-validated AUC at lambda.min
    lambda_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
    auc_cv <- cv_fit$cvm[lambda_idx]

    if (verbose) {
        cat(sprintf("Cross-validated AUC: %.3f\n", auc_cv))
    }

    ## ---- pathway importance from glmnet coefficients -------------------------
    coefs <- as.numeric(
        stats::coef(cv_fit, s = "lambda.min")[-1]  # drop intercept
    )
    names(coefs) <- pretrained@pathway_names
    pathway_importance <- sort(abs(coefs), decreasing = TRUE)

    if (verbose) {
        n_nonzero <- sum(coefs != 0)
        cat(sprintf("Non-zero pathway coefficients: %d / %d\n",
                    n_nonzero, length(coefs)))
        top5 <- head(pathway_importance, 5)
        cat("Top 5 pathways:\n")
        for (i in seq_along(top5)) {
            cat(sprintf("  %d. %s (|coef| = %.4f)\n",
                        i, names(top5)[i], top5[i]))
        }
    }

    ## ---- return results ------------------------------------------------------
    list(
        predictions         = preds,
        labels              = y_factor,
        auc                 = auc_cv,
        pathway_activations = H,
        pathway_importance  = pathway_importance,
        glmnet_fit          = cv_fit,
        gene_coverage       = aligned$n_matched / nrow(pretrained@gene_info),
        lambda              = cv_fit$lambda.min,
        class_levels        = levels_y
    )
}


# ── Utility: extract pretrained init for Mode 2 ──────────────────

#' Extract Pretrained Weights Aligned to an Architecture
#'
#' Aligns pretrained hidden-layer weights to the gene order defined by
#' a \code{VNNArchitecture}'s first-layer mask. Used internally by
#' \code{\link{trainVNN}} when the \code{pretrained} argument is set
#' (Mode 2: fine-tuned transfer learning).
#'
#' @param pretrained A \code{PretrainedVNN} object.
#' @param architecture A \code{VNNArchitecture} object. The first-layer
#'   mask must have rownames matching the pretrained gene IDs.
#' @param gene_id_type One of \code{"symbol"}, \code{"ensembl"}, or
#'   \code{"auto"}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{W_init}}{Dense matrix \code{[n_arch_genes x n_pathways]}
#'       of aligned pretrained weights. Genes not in the pretrained model
#'       receive zero initialization.}
#'     \item{\code{b_init}}{Numeric vector of pretrained biases.}
#'     \item{\code{n_matched}}{Number of architecture genes found in
#'       pretrained model.}
#'   }
#'
#' @keywords internal
.alignPretrainedToArchitecture <- function(pretrained, architecture,
                                           gene_id_type = "auto") {

    masks <- layerMasks(architecture)
    arch_mask <- masks[[1]]
    arch_genes <- rownames(arch_mask)

    if (is.null(arch_genes)) {
        stop("Architecture first-layer mask must have rownames (gene IDs).",
             call. = FALSE)
    }

    ## Detect ID type
    if (gene_id_type == "auto") {
        gene_id_type <- if (any(grepl("^ENSG", head(arch_genes, 20)))) {
            "ensembl"
        } else {
            "symbol"
        }
    }

    pt_col <- switch(gene_id_type,
        symbol  = "symbol",
        ensembl = "ensembl"
    )
    pt_ids <- pretrained@gene_info[[pt_col]]

    ## Strip Ensembl versions
    arch_genes_clean <- arch_genes
    if (gene_id_type == "ensembl") {
        arch_genes_clean <- sub("\\.\\d+$", "", arch_genes_clean)
    }

    ## Build aligned weight matrix (zeros for unmatched genes)
    n_arch_genes <- length(arch_genes)
    n_pathways   <- ncol(pretrained@W_hidden)

    W_init <- matrix(0, nrow = n_arch_genes, ncol = n_pathways,
                     dimnames = list(arch_genes, pretrained@pathway_names))

    pt_lookup <- setNames(seq_along(pt_ids), pt_ids)
    matched <- intersect(arch_genes_clean, pt_ids)

    for (g in matched) {
        arch_row <- which(arch_genes_clean == g)[1]
        pt_row   <- pt_lookup[g]
        W_init[arch_row, ] <- pretrained@W_hidden[pt_row, ]
    }

    message(sprintf(
        "Pretrained init: %d / %d architecture genes matched (%.1f%%)",
        length(matched), n_arch_genes,
        length(matched) / n_arch_genes * 100
    ))

    list(
        W_init    = W_init,
        b_init    = pretrained@b_hidden,
        n_matched = length(matched)
    )
}
