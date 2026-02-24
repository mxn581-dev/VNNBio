#' @title VNNBio: Visible Neural Networks for Biology
#'
#' @description Implements Visible Neural Networks (VNNs) whose architecture is
#'   constrained by prior biological knowledge such as gene-pathway mappings.
#'   Network connectivity mirrors known biology, so each hidden node corresponds
#'   to an interpretable biological concept.
#'
#' @import methods
#' @importFrom Matrix sparseMatrix nnzero
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'   colData
#' @importFrom S4Vectors DataFrame
#' @importFrom stats predict setNames
#' @importFrom ggplot2 ggplot aes geom_col geom_hline coord_flip labs
#'   theme_minimal theme element_text scale_fill_gradient2 facet_wrap
#'   geom_segment .data
#' @importFrom utils head tail
#'
#' @aliases VNNBio-package
#' @keywords internal
"_PACKAGE"
