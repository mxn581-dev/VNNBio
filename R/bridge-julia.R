# =============================================================================
# bridge-julia.R -- R-to-Julia Bridge for VNN Training
# =============================================================================

# =============================================================================
# Module-level state (package environment)
# =============================================================================

.vnn_env <- new.env(parent = emptyenv())
.vnn_env$julia_ready   <- FALSE
.vnn_env$vnn_module     <- NULL


# =============================================================================
# Julia initialization (called from .onLoad in zzz.R)
# =============================================================================

#' Initialize the Julia Backend
#'
#' Activates the Julia project bundled in \code{inst/julia/}, installs
#' dependencies if needed, and loads the VNNCore module.
#'
#' @param julia_path Path to the Julia binary. If NULL, uses the system
#'   default found by JuliaConnectoR.
#' @param force Logical. Re-initialize even if already loaded? Default FALSE.
#'
#' @return Invisible TRUE on success.
#' @export
initJuliaBackend <- function(julia_path = NULL, force = FALSE) {

    if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
        stop("Package 'JuliaConnectoR' is required for model training.\n",
             "Install with: install.packages('JuliaConnectoR')",
             call. = FALSE)
    }

    if (.vnn_env$julia_ready && !force) {
        message("Julia backend already initialized. Use force=TRUE to reload.")
        return(invisible(TRUE))
    }

    ## Locate the Julia project shipped with this package
    julia_proj <- system.file("julia", package = "VNNBio", mustWork = TRUE)

    ## Tell Julia to use our bundled project
    JuliaConnectoR::juliaEval(paste0('
        import Pkg
        Pkg.activate("', julia_proj, '")
        Pkg.instantiate()
    '))

    ## Load VNNCore module
    JuliaConnectoR::juliaEval(paste0(
        'include("', file.path(julia_proj, "src", "VNNCore.jl"), '")'
    ))
    JuliaConnectoR::juliaEval("using .VNNCore")

    ## Julia 1.12+ requires invokelatest for functions defined via
    ## include() in a different world age. Define top-level wrappers
    ## so that JuliaConnectoR calls go through invokelatest.
    JuliaConnectoR::juliaEval('
        _vnn_receive_sparse_mask(args...) = Base.invokelatest(VNNCore.receive_sparse_mask, args...)
        _vnn_receive_data(args...)        = Base.invokelatest(VNNCore.receive_data, args...)
        _vnn_train_vnn(args...)           = Base.invokelatest(VNNCore.train_vnn, args...)
        _vnn_get_masked_weights(args...)  = Base.invokelatest(VNNCore.get_masked_weights, args...)
        _vnn_predict_vnn(args...)         = Base.invokelatest(VNNCore.predict_vnn, args...)
    ')

    .vnn_env$julia_ready <- TRUE
    message("VNNBio Julia backend initialized.")
    invisible(TRUE)
}


# =============================================================================
# Data serialization: R sparse matrix -> Julia SparseMatrixCSC
# =============================================================================

#' Transfer a Sparse Mask to Julia
#'
#' Converts an R \code{dgCMatrix} into Julia's native
#' \code{SparseMatrixCSC{Float32, Int64}} format by passing the raw CSC
#' triplet components. This avoids materializing a dense matrix.
#'
#' @param mask A \code{dgCMatrix} from the \code{Matrix} package.
#' @param julia_varname Character. Name to assign in Julia's Main scope.
#'
#' @return Invisible Julia reference.
#' @keywords internal
.transferSparseMask <- function(mask, julia_varname) {

    stopifnot(inherits(mask, "dgCMatrix"))

    ## dgCMatrix stores: @i (0-based row indices), @p (column pointers), @x
    ## Julia SparseMatrixCSC uses 1-based indices
    row_idx <- as.integer(mask@i + 1L)       # 0-based -> 1-based
    col_ptr <- as.integer(mask@p + 1L)       # 0-based -> 1-based
    values  <- as.numeric(mask@x)
    dims    <- dim(mask)

    ## Ship the three vectors + dims to Julia, construct sparse matrix there
    JuliaConnectoR::juliaCall("_vnn_receive_sparse_mask",
        row_idx, col_ptr, values,
        as.integer(dims[1]), as.integer(dims[2]),
        julia_varname
    )
}


# =============================================================================
# Training dispatch
# =============================================================================

#' Train a VNN Model via Julia
#'
#' Serializes the expression matrix, label vector, and architecture masks to
#' Julia, then calls the Lux.jl training loop.
#'
#' @param se A \code{SummarizedExperiment} with the expression assay.
#' @param architecture A \code{VNNArchitecture} object.
#' @param label_col Character. Column in \code{colData(se)} to use as the
#'   response variable.
#' @param assay_name Character. Name of the assay to use. Default "counts" or
#'   first available.
#' @param task One of "classification", "regression".
#' @param epochs Integer. Number of training epochs. Default 100.
#' @param learning_rate Numeric. Adam optimizer LR. Default 1e-3.
#' @param batch_size Integer. Mini-batch size. Default 64.
#' @param val_fraction Numeric in (0,1). Fraction held out for validation.
#' @param l1_lambda Numeric. L1 penalty on masked weights for sparsity.
#' @param patience Integer. Stop training if validation loss does not improve
#'   for this many consecutive epochs. Set to 0 to disable early stopping.
#'   Default 0 (disabled).
#' @param min_delta Numeric. Minimum improvement in validation loss to be
#'   considered progress. Default 1e-4.
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print per-epoch progress? Default TRUE.
#'
#' @return A \code{VNNModel} object.
#' @export
trainVNN <- function(se,
                     architecture,
                     label_col,
                     assay_name = NULL,
                     task = c("classification", "regression"),
                     epochs = 100L,
                     learning_rate = 1e-3,
                     batch_size = 64L,
                     val_fraction = 0.2,
                     l1_lambda = 1e-4,
                     patience = 0L,
                     min_delta = 1e-4,
                     seed = 42L,
                     verbose = TRUE) {

    task <- match.arg(task)

    ## ---- ensure Julia is ready ----------------------------------------------
    if (!.vnn_env$julia_ready) initJuliaBackend()

    ## ---- extract expression matrix (genes x samples) -> (samples x genes) ---
    if (is.null(assay_name)) {
        assay_name <- SummarizedExperiment::assayNames(se)[1]
    }
    X <- t(as.matrix(SummarizedExperiment::assay(se, assay_name)))
    ## Now X is [n_samples x n_genes]

    ## ---- extract labels -----------------------------------------------------
    y <- SummarizedExperiment::colData(se)[[label_col]]
    if (task == "classification") {
        y <- as.integer(as.factor(y)) - 1L  # 0-indexed for binary cross-entropy
    }
    y <- as.numeric(y)

    ## ---- transfer masks to Julia --------------------------------------------
    masks <- layerMasks(architecture)
    mask_names <- character(length(masks))
    for (i in seq_along(masks)) {
        varname <- paste0("mask_layer_", i)
        .transferSparseMask(masks[[i]], varname)
        mask_names[i] <- varname
    }

    ## ---- transfer data to Julia ---------------------------------------------
    JuliaConnectoR::juliaCall("_vnn_receive_data",
        X, y,
        as.integer(nOutput(architecture)),
        activationFunction(architecture)
    )

    ## ---- run training loop --------------------------------------------------
    history_raw <- JuliaConnectoR::juliaCall("_vnn_train_vnn",
        mask_names,
        as.integer(epochs),
        learning_rate,
        as.integer(batch_size),
        val_fraction,
        l1_lambda,
        as.integer(patience),
        min_delta,
        as.integer(seed),
        verbose
    )

    ## ---- pull back training history -----------------------------------------
    ## Early stopping may return fewer epochs than requested
    train_losses <- as.numeric(history_raw$train_losses)
    val_losses <- as.numeric(history_raw$val_losses)
    actual_epochs <- length(train_losses)

    history_df <- data.frame(
        epoch      = seq_len(actual_epochs),
        train_loss = train_losses,
        val_loss   = val_losses
    )

    ## ---- extract importance scores ------------------------------------------
    importance <- .extractImportanceScores(architecture,
                                           history_raw$model_ref)

    ## ---- assemble VNNModel --------------------------------------------------
    model <- new("VNNModel",
        architecture      = architecture,
        julia_model_ref   = history_raw$model_ref,
        training_history  = history_df,
        importance_scores = importance,
        input_se          = se,
        task              = task
    )

    model
}


# =============================================================================
# Importance extraction
# =============================================================================

#' Extract Node-Level Importance Scores from Julia
#'
#' After training, queries the Julia model for each layer's weight matrix,
#' applies the mask, and computes per-node importance as the column-wise sum
#' of absolute (masked) weights.
#'
#' @param architecture A \code{VNNArchitecture} object (for node names).
#' @return Named list of named numeric vectors.
#' @keywords internal
.extractImportanceScores <- function(architecture, model_ref = "") {

    scores <- list()
    masks <- layerMasks(architecture)
    lnames <- layerNames(architecture)

    for (i in seq_along(masks)) {
        ## Ask Julia for the masked weight matrix of layer i
        raw <- JuliaConnectoR::juliaCall("_vnn_get_masked_weights", as.integer(i),
                          model_ref)
        W <- as.matrix(raw)  # [in_features x out_features]

        ## Per-output-node importance = sum of |w_ij| over input dimension
        node_importance <- colSums(abs(W))

        ## Attach pathway/node names from the mask column names
        mask_i <- masks[[i]]
        if (!is.null(colnames(mask_i))) {
            names(node_importance) <- colnames(mask_i)
        }

        layer_label <- if (length(lnames) >= i) lnames[i] else paste0("layer_", i)
        scores[[layer_label]] <- sort(node_importance, decreasing = TRUE)
    }

    scores
}
