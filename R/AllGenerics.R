# =============================================================================
# AllGenerics.R -- S4 Generic Function Definitions
#
# Bioconductor convention: ALL generics in one file, alphabetically ordered.
# Accessors follow the pattern: slot name as getter, `slot<-` as setter.
# =============================================================================


# -- Accessors: GenePathwayMap ------------------------------------

#' @rdname GenePathwayMap-class
#' @export
setGeneric("mappingTable", function(object) standardGeneric("mappingTable"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("maskMatrix", function(object) standardGeneric("maskMatrix"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("geneIndex", function(object) standardGeneric("geneIndex"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("pathwayIndex", function(object) standardGeneric("pathwayIndex"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("maskSource", function(object) standardGeneric("maskSource"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("maskDensity", function(object) standardGeneric("maskDensity"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("nGenes", function(object) standardGeneric("nGenes"))

#' @rdname GenePathwayMap-class
#' @export
setGeneric("nPathways", function(object) standardGeneric("nPathways"))


# -- Accessors: VNNArchitecture -----------------------------------

#' @rdname VNNArchitecture-class
#' @export
setGeneric("layerMasks", function(object) standardGeneric("layerMasks"))

#' @rdname VNNArchitecture-class
#' @export
setGeneric("layerNames", function(object) standardGeneric("layerNames"))

#' @rdname VNNArchitecture-class
#' @export
setGeneric("activationFunction", function(object)
    standardGeneric("activationFunction"))

#' @rdname VNNArchitecture-class
#' @export
setGeneric("nOutput", function(object) standardGeneric("nOutput"))

#' @rdname VNNArchitecture-class
#' @export
setGeneric("nLayers", function(object) standardGeneric("nLayers"))

#' @rdname VNNArchitecture-class
#' @export
setGeneric("validateChain", function(object) standardGeneric("validateChain"))


# -- Accessors: VNNModel ------------------------------------------

#' @rdname VNNModel-class
#' @export
setGeneric("architecture", function(object) standardGeneric("architecture"))

#' @rdname VNNModel-class
#' @export
setGeneric("trainingHistory", function(object)
    standardGeneric("trainingHistory"))

#' @rdname VNNModel-class
#' @export
setGeneric("importanceScores", function(object)
    standardGeneric("importanceScores"))

#' @rdname VNNModel-class
#' @export
setGeneric("taskType", function(object) standardGeneric("taskType"))

#' @rdname VNNModel-class
#' @export
setGeneric("inputSE", function(object) standardGeneric("inputSE"))


# -- Core methods -------------------------------------------------

#' Construct a VNNArchitecture from a GenePathwayMap
#'
#' @param object A \code{GenePathwayMap} object.
#' @param ... Additional arguments passed to constructors.
#'
#' @return A \code{VNNArchitecture} object.
#' @export
setGeneric("buildArchitecture", function(object, ...)
    standardGeneric("buildArchitecture"))

#' Extract pathway importance as a tidy data.frame
#'
#' @param object A \code{VNNModel} object.
#' @param ... Additional arguments.
#'
#' @return A \code{data.frame} with columns \code{layer}, \code{pathway},
#'   \code{importance}, and \code{rank}.
#' @export
setGeneric("importanceTable", function(object, ...)
    standardGeneric("importanceTable"))
