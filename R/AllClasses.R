# =============================================================================
# AllClasses.R -- S4 Class Definitions
#
# Bioconductor convention: ALL class definitions in one file.
# Each class has: full roxygen docs, slots, prototype, validity.
# =============================================================================

# -- GenePathwayMap -----------------------------------------------

#' GenePathwayMap: Biological Prior Knowledge as a Sparse Matrix
#'
#' Stores the mapping from genes to biological pathways and the derived sparse
#' adjacency mask used to constrain neural network connectivity. This is the
#' central object that encodes prior biological knowledge into the VNN
#' architecture.
#'
#' @slot mapping A \code{data.frame} with at minimum columns \code{gene_id}
#'   and \code{pathway_id}. May optionally include \code{gene_symbol},
#'   \code{pathway_name}, and a numeric \code{weight} column.
#' @slot mask A sparse binary matrix (\code{dgCMatrix}) of dimensions
#'   \code{[n_genes x n_pathways]}. Entry \code{(i, j) = 1} if gene
#'   \code{i} belongs to pathway \code{j}.
#' @slot gene_index Named character vector mapping gene IDs to their row
#'   position in the mask matrix.
#' @slot pathway_index Named character vector mapping pathway IDs to their
#'   column position in the mask matrix.
#' @slot source Character scalar indicating the ontology source
#'   (e.g., \code{"KEGG"}, \code{"GO_BP"}, \code{"MSigDB_H"}).
#'
#' @section Constructors:
#' Use \code{\link{buildGenePathwayMap}} or \code{\link{buildMapFromMSigDB}}
#' rather than calling \code{new("GenePathwayMap", ...)} directly.
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{mappingTable(x)}}{Returns the gene-pathway data.frame.}
#'   \item{\code{maskMatrix(x)}}{Returns the sparse mask (dgCMatrix).}
#'   \item{\code{geneIndex(x)}}{Returns the ordered gene ID vector.}
#'   \item{\code{pathwayIndex(x)}}{Returns the ordered pathway ID vector.}
#'   \item{\code{maskSource(x)}}{Returns the ontology source label.}
#'   \item{\code{nGenes(x)}}{Returns the number of genes (mask rows).}
#'   \item{\code{nPathways(x)}}{Returns the number of pathways (mask cols).}
#'   \item{\code{maskDensity(x)}}{Returns the fraction of nonzero entries.}
#' }
#'
#' @seealso \code{\link{buildGenePathwayMap}}, \code{\link{buildMapFromMSigDB}}
#'
#' @param object A \code{GenePathwayMap} object.
#' @param x A \code{GenePathwayMap} object (for S4 methods).
#' @param ... Additional arguments (unused).
#'
#' @name GenePathwayMap-class
#' @rdname GenePathwayMap-class
#' @exportClass GenePathwayMap
#' @examples
#' \donttest{
#' mapping_df <- data.frame(
#'     gene_id    = c("G1", "G1", "G2", "G3", "G4", "G5"),
#'     pathway_id = c("P_A", "P_B", "P_A", "P_A", "P_B", "P_B")
#' )
#' gpm <- buildGenePathwayMap(mapping_df, min_pathway_size = 2)
#' maskDensity(gpm)
#' nGenes(gpm)
#' }
setClass("GenePathwayMap",
    slots = list(
        mapping       = "data.frame",
        mask          = "dgCMatrix",
        gene_index    = "character",
        pathway_index = "character",
        source        = "character"
    ),
    prototype = list(
        mapping       = data.frame(gene_id = character(0),
                                   pathway_id = character(0)),
        mask          = new("dgCMatrix"),
        gene_index    = character(0),
        pathway_index = character(0),
        source        = "custom"
    )
)

setValidity("GenePathwayMap", function(object) {
    errors <- character()
    req_cols <- c("gene_id", "pathway_id")
    if (!all(req_cols %in% colnames(object@mapping))) {
        errors <- c(errors,
            paste0("'mapping' must contain columns: ",
                   paste(req_cols, collapse = ", ")))
    }
    if (nrow(object@mask) > 0) {
        if (nrow(object@mask) != length(object@gene_index))
            errors <- c(errors,
                "'mask' row count must equal length of 'gene_index'")
        if (ncol(object@mask) != length(object@pathway_index))
            errors <- c(errors,
                "'mask' column count must equal length of 'pathway_index'")
    }
    if (length(object@source) != 1 || is.na(object@source))
        errors <- c(errors, "'source' must be a single non-NA string")

    if (length(errors) == 0) TRUE else errors
})


# -- VNNArchitecture ----------------------------------------------

#' VNNArchitecture: Layer-wise Network Topology
#'
#' Encodes the full Visible Neural Network topology: the ordered list of
#' sparse masks from input genes through intermediate pathway layers to the
#' output layer. Supports hierarchical ontologies (e.g., GO terms organized
#' in a DAG) represented as a list of layer masks.
#'
#' @slot layer_masks A named list of \code{dgCMatrix} objects. Each matrix
#'   defines connectivity from one layer to the next. Names should follow
#'   the pattern \code{"input_to_pathway"}, \code{"pathway_to_system"}, etc.
#' @slot layer_names Character vector of human-readable layer names.
#' @slot activation Character scalar: activation function applied at hidden
#'   nodes. One of \code{"tanh"}, \code{"relu"}, \code{"sigmoid"},
#'   \code{"gelu"}, \code{"swish"}.
#' @slot n_output Integer: number of output nodes (1 for binary
#'   classification, k for multi-class).
#'
#' @section Validity:
#' For multi-layer architectures, the number of columns in mask \code{L}
#' must equal the number of rows in mask \code{L+1}. Use
#' \code{validateChain()} to check this explicitly.
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{layerMasks(x)}}{Returns the list of sparse mask matrices.}
#'   \item{\code{layerNames(x)}}{Returns the character vector of layer names.}
#'   \item{\code{activationFunction(x)}}{Returns the activation name.}
#'   \item{\code{nOutput(x)}}{Returns the number of output nodes.}
#'   \item{\code{nLayers(x)}}{Returns the number of masked layers.}
#'   \item{\code{validateChain(x)}}{Checks dimensional compatibility.}
#' }
#'
#' @seealso \code{\link{buildArchitecture}}
#'
#' @param object A \code{VNNArchitecture} object.
#'
#' @name VNNArchitecture-class
#' @rdname VNNArchitecture-class
#' @exportClass VNNArchitecture
setClass("VNNArchitecture",
    slots = list(
        layer_masks  = "list",
        layer_names  = "character",
        activation   = "character",
        n_output     = "integer"
    ),
    prototype = list(
        layer_masks  = list(),
        layer_names  = character(0),
        activation   = "tanh",
        n_output     = 1L
    )
)

setValidity("VNNArchitecture", function(object) {
    errors <- character()
    valid_act <- c("tanh", "relu", "sigmoid", "gelu", "swish")
    if (length(object@activation) != 1 ||
        !object@activation %in% valid_act) {
        errors <- c(errors,
            paste0("'activation' must be one of: ",
                   paste(valid_act, collapse = ", ")))
    }
    if (length(object@n_output) != 1 || object@n_output < 1L)
        errors <- c(errors, "'n_output' must be a positive integer")

    ## Check layer mask types
    for (i in seq_along(object@layer_masks)) {
        if (!is(object@layer_masks[[i]], "dgCMatrix") &&
            !is(object@layer_masks[[i]], "dMatrix"))
            errors <- c(errors,
                paste0("layer_masks[[", i, "]] must be a sparse Matrix"))
    }

    ## Check dimensional chaining
    masks <- object@layer_masks
    if (length(masks) > 1) {
        for (i in seq_len(length(masks) - 1)) {
            if (ncol(masks[[i]]) != nrow(masks[[i + 1]])) {
                errors <- c(errors,
                    sprintf(paste0("Dimension mismatch: ncol(layer %d) = %d ",
                                   "but nrow(layer %d) = %d"),
                            i, ncol(masks[[i]]),
                            i + 1, nrow(masks[[i + 1]])))
            }
        }
    }

    if (length(errors) == 0) TRUE else errors
})


# -- VNNModel -----------------------------------------------------

#' VNNModel: Trained Visible Neural Network
#'
#' Stores the trained VNN model including Julia-side parameter references,
#' architecture specification, training history, and extracted
#' interpretability scores.
#'
#' @slot architecture A \code{\link{VNNArchitecture-class}} object.
#' @slot julia_model_ref An external reference to the Julia-side trained
#'   model parameters (managed by JuliaConnectoR). Users should not
#'   interact with this directly.
#' @slot training_history A \code{data.frame} with columns \code{epoch},
#'   \code{train_loss}, \code{val_loss}.
#' @slot importance_scores A named list of named numeric vectors, one per
#'   biological layer, containing node-level importance scores (sum of
#'   absolute masked weights).
#' @slot input_se The original \code{SummarizedExperiment} used for training,
#'   stored for provenance and feature alignment at prediction time.
#' @slot task Character: one of \code{"classification"}, \code{"regression"},
#'   \code{"regression"}.
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{architecture(x)}}{Returns the VNNArchitecture.}
#'   \item{\code{trainingHistory(x)}}{Returns training loss data.frame.}
#'   \item{\code{importanceScores(x)}}{Returns the named list of scores.}
#'   \item{\code{importanceTable(x)}}{Returns a tidy data.frame of scores.}
#'   \item{\code{taskType(x)}}{Returns the task type string.}
#'   \item{\code{inputSE(x)}}{Returns the input SummarizedExperiment.}
#'   \item{\code{predict(x, newdata)}}{Predict on new data.}
#' }
#'
#' @name VNNModel-class
#' @rdname VNNModel-class
#' @exportClass VNNModel
setClass("VNNModel",
    slots = list(
        architecture      = "VNNArchitecture",
        julia_model_ref   = "ANY",
        training_history  = "data.frame",
        importance_scores = "list",
        input_se          = "ANY",
        task              = "character"
    ),
    prototype = list(
        architecture      = new("VNNArchitecture"),
        julia_model_ref   = NULL,
        training_history  = data.frame(epoch = integer(0),
                                       train_loss = numeric(0),
                                       val_loss = numeric(0)),
        importance_scores = list(),
        input_se          = NULL,
        task              = "classification"
    )
)

setValidity("VNNModel", function(object) {
    errors <- character()
    valid_tasks <- c("classification", "regression")
    if (length(object@task) != 1 || !object@task %in% valid_tasks)
        errors <- c(errors,
            paste0("'task' must be one of: ",
                   paste(valid_tasks, collapse = ", ")))

    hist_cols <- c("epoch", "train_loss", "val_loss")
    if (nrow(object@training_history) > 0 &&
        !all(hist_cols %in% colnames(object@training_history)))
        errors <- c(errors,
            paste0("'training_history' must contain columns: ",
                   paste(hist_cols, collapse = ", ")))

    if (length(errors) == 0) TRUE else errors
})
