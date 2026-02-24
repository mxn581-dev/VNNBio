# =============================================================================
# AllMethods.R -- S4 Method Implementations
#
# Bioconductor convention: ALL method implementations in one file.
# Organized by class, then alphabetically within each class.
# =============================================================================


# =========================================================================== #
#                          GenePathwayMap methods                              #
# =========================================================================== #

#' @rdname GenePathwayMap-class
#' @aliases show,GenePathwayMap-method
setMethod("show", "GenePathwayMap", function(object) {
    n_genes <- nrow(object@mask)
    n_pw    <- ncol(object@mask)
    nnz     <- nnzero(object@mask)
    density <- if (n_genes * n_pw > 0)
        round(100 * nnz / (n_genes * n_pw), 2) else 0

    cat("GenePathwayMap object\n")
    cat("  Source:    ", object@source, "\n")
    cat("  Genes:     ", n_genes, "\n")
    cat("  Pathways:  ", n_pw, "\n")
    cat("  Nonzeros:  ", nnz, " (", density, "% density)\n", sep = "")
    cat("  Mapping:   ", nrow(object@mapping), " gene-pathway pairs\n")
    if (n_pw > 0) {
        n_show <- min(5, n_pw)
        cat("  Top pathways: ",
            paste(head(colnames(object@mask), n_show), collapse = ", "),
            if (n_pw > n_show) ", ..." else "", "\n")
    }
})

#' @rdname GenePathwayMap-class
#' @aliases mappingTable,GenePathwayMap-method
setMethod("mappingTable", "GenePathwayMap", function(object) object@mapping)

#' @rdname GenePathwayMap-class
#' @aliases maskMatrix,GenePathwayMap-method
setMethod("maskMatrix", "GenePathwayMap", function(object) object@mask)

#' @rdname GenePathwayMap-class
#' @aliases geneIndex,GenePathwayMap-method
setMethod("geneIndex", "GenePathwayMap", function(object) object@gene_index)

#' @rdname GenePathwayMap-class
#' @aliases pathwayIndex,GenePathwayMap-method
setMethod("pathwayIndex", "GenePathwayMap", function(object)
    object@pathway_index)

#' @rdname GenePathwayMap-class
#' @aliases maskSource,GenePathwayMap-method
setMethod("maskSource", "GenePathwayMap", function(object) object@source)

#' @rdname GenePathwayMap-class
#' @aliases nGenes,GenePathwayMap-method
setMethod("nGenes", "GenePathwayMap", function(object) nrow(object@mask))

#' @rdname GenePathwayMap-class
#' @aliases nPathways,GenePathwayMap-method
setMethod("nPathways", "GenePathwayMap", function(object) ncol(object@mask))

#' @rdname GenePathwayMap-class
#' @aliases maskDensity,GenePathwayMap-method
setMethod("maskDensity", "GenePathwayMap", function(object) {
    total <- nrow(object@mask) * ncol(object@mask)
    if (total == 0) return(0)
    nnzero(object@mask) / total
})

#' @rdname GenePathwayMap-class
#' @aliases dim,GenePathwayMap-method
#' @export
setMethod("dim", "GenePathwayMap", function(x) dim(x@mask))

#' @rdname GenePathwayMap-class
#' @aliases length,GenePathwayMap-method
#' @export
setMethod("length", "GenePathwayMap", function(x) nrow(x@mapping))


# =========================================================================== #
#                         VNNArchitecture methods                              #
# =========================================================================== #

#' @rdname VNNArchitecture-class
#' @aliases show,VNNArchitecture-method
setMethod("show", "VNNArchitecture", function(object) {
    n <- length(object@layer_masks)
    cat("VNNArchitecture object\n")
    cat("  Layers:     ", n, "\n")
    cat("  Activation: ", object@activation, "\n")
    cat("  Output dim: ", object@n_output, "\n")
    if (n > 0) {
        cat("  Layer topology:\n")
        for (i in seq_len(n)) {
            mask <- object@layer_masks[[i]]
            lname <- if (i <= length(object@layer_names))
                object@layer_names[i] else paste0("layer_", i)
            cat(sprintf("    [%d] %s: %d -> %d (%.1f%% sparse)\n",
                i, lname, nrow(mask), ncol(mask),
                100 * (1 - nnzero(mask) / (nrow(mask) * ncol(mask)))))
        }
        last_mask <- object@layer_masks[[n]]
        cat(sprintf("    [%d] output: %d -> %d (dense)\n",
            n + 1, ncol(last_mask), object@n_output))
    }
})

#' @rdname VNNArchitecture-class
#' @aliases layerMasks,VNNArchitecture-method
setMethod("layerMasks", "VNNArchitecture", function(object)
    object@layer_masks)

#' @rdname VNNArchitecture-class
#' @aliases layerNames,VNNArchitecture-method
setMethod("layerNames", "VNNArchitecture", function(object)
    object@layer_names)

#' @rdname VNNArchitecture-class
#' @aliases activationFunction,VNNArchitecture-method
setMethod("activationFunction", "VNNArchitecture", function(object)
    object@activation)

#' @rdname VNNArchitecture-class
#' @aliases nOutput,VNNArchitecture-method
setMethod("nOutput", "VNNArchitecture", function(object) object@n_output)

#' @rdname VNNArchitecture-class
#' @aliases nLayers,VNNArchitecture-method
setMethod("nLayers", "VNNArchitecture", function(object)
    length(object@layer_masks))

#' @rdname VNNArchitecture-class
#' @aliases validateChain,VNNArchitecture-method
setMethod("validateChain", "VNNArchitecture", function(object) {
    masks <- object@layer_masks
    if (length(masks) <= 1) return(TRUE)
    for (i in seq_len(length(masks) - 1)) {
        if (ncol(masks[[i]]) != nrow(masks[[i + 1]])) {
            stop(sprintf(
                paste0("Layer dimension mismatch: layer %d outputs %d ",
                       "nodes but layer %d expects %d inputs"),
                i, ncol(masks[[i]]), i + 1, nrow(masks[[i + 1]])))
        }
    }
    TRUE
})


# -- buildArchitecture convenience method -------------------------

#' @rdname GenePathwayMap-class
#' @aliases buildArchitecture,GenePathwayMap-method
#'
#' @param activation Character. Activation function name. Default "tanh".
#' @param n_output Integer. Number of output nodes. Default 1L.
#'
#' @examples
#' mapping_df <- data.frame(
#'     gene_id    = c("G1","G1","G2","G3","G4","G5"),
#'     pathway_id = c("PA","PB","PA","PA","PB","PB")
#' )
#' gpm <- buildGenePathwayMap(mapping_df, min_pathway_size = 2)
#' arch <- buildArchitecture(gpm)
setMethod("buildArchitecture", "GenePathwayMap",
    function(object, activation = "tanh", n_output = 1L, ...) {
        new("VNNArchitecture",
            layer_masks = list(input_to_pathway = object@mask),
            layer_names = paste0(object@source, " pathways"),
            activation  = activation,
            n_output    = as.integer(n_output)
        )
    }
)


# =========================================================================== #
#                            VNNModel methods                                  #
# =========================================================================== #

#' @rdname VNNModel-class
#' @aliases show,VNNModel-method
setMethod("show", "VNNModel", function(object) {
    n_epochs <- nrow(object@training_history)
    cat("VNNModel object\n")
    cat("  Task:        ", object@task, "\n")
    cat("  Epochs:      ", n_epochs, "\n")
    if (n_epochs > 0) {
        cat("  Final loss:   train=",
            round(tail(object@training_history$train_loss, 1), 4),
            "  val=",
            round(tail(object@training_history$val_loss, 1), 4), "\n")
    }
    cat("  Architecture: ")
    show(object@architecture)
    if (length(object@importance_scores) > 0) {
        cat("  Importance scores available for ",
            length(object@importance_scores), " layer(s)\n")
        for (nm in names(object@importance_scores)) {
            s <- object@importance_scores[[nm]]
            n_show <- min(5, length(s))
            cat("    ", nm, ": ",
                paste(paste0(names(s)[seq_len(n_show)], " (",
                             round(s[seq_len(n_show)], 3), ")"),
                      collapse = ", "),
                if (length(s) > n_show) ", ..." else "", "\n", sep = "")
        }
    }
})

#' @rdname VNNModel-class
#' @aliases architecture,VNNModel-method
setMethod("architecture", "VNNModel", function(object) object@architecture)

#' @rdname VNNModel-class
#' @aliases trainingHistory,VNNModel-method
setMethod("trainingHistory", "VNNModel", function(object)
    object@training_history)

#' @rdname VNNModel-class
#' @aliases importanceScores,VNNModel-method
setMethod("importanceScores", "VNNModel", function(object)
    object@importance_scores)

#' @rdname VNNModel-class
#' @aliases taskType,VNNModel-method
setMethod("taskType", "VNNModel", function(object) object@task)

#' @rdname VNNModel-class
#' @aliases inputSE,VNNModel-method
setMethod("inputSE", "VNNModel", function(object) object@input_se)

#' @rdname VNNModel-class
#' @aliases importanceTable,VNNModel-method
setMethod("importanceTable", "VNNModel", function(object, ...) {
    scores <- object@importance_scores
    if (length(scores) == 0)
        return(data.frame(layer = character(0), pathway = character(0),
                          importance = numeric(0), rank = integer(0)))
    do.call(rbind, lapply(names(scores), function(layer) {
        s <- scores[[layer]]
        data.frame(
            layer      = layer,
            pathway    = names(s),
            importance = unname(s),
            rank       = seq_along(s),
            stringsAsFactors = FALSE
        )
    }))
})


# -- predict method -----------------------------------------------

#' Predict Using a Trained VNN Model
#'
#' @param object A trained \code{VNNModel}.
#' @param newdata A \code{SummarizedExperiment} with the same features
#'   (genes) as the training data, or a numeric matrix in
#'   \code{[samples x genes]} orientation (rows = samples).
#' @param assay_name Character. Which assay to use from the SE.
#'   Default: first available.
#' @param ... Additional arguments (unused).
#'
#' @return A numeric vector of predictions. For classification, values
#'   are probabilities between 0 and 1.
#'
#' @rdname VNNModel-class
#' @aliases predict,VNNModel-method
#' @exportMethod predict
setMethod("predict", "VNNModel",
    function(object, newdata, assay_name = NULL, ...) {
        model_ref <- object@julia_model_ref
        if (is.null(model_ref))
            stop("Model has no trained parameters. Run trainVNN() first.")

        if (is(newdata, "SummarizedExperiment")) {
            if (is.null(assay_name))
                assay_name <- SummarizedExperiment::assayNames(newdata)[1]
            ## SE assays are [genes x samples]; Julia expects [samples x genes]
            X <- t(as.matrix(SummarizedExperiment::assay(newdata,
                                                          assay_name)))
        } else if (is.matrix(newdata)) {
            ## Require explicit [samples x genes] -- never guess orientation
            X <- newdata
        } else {
            stop("'newdata' must be a SummarizedExperiment or a numeric ",
                 "matrix with rows = samples and columns = genes")
        }

        preds <- JuliaConnectoR::juliaCall("_vnn_predict_vnn", X, model_ref)
        as.numeric(preds)
    }
)
