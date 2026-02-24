# =============================================================================
# mapping-utils.R -- Gene-Pathway Mapping Utilities
# =============================================================================

# =============================================================================
# Core: data.frame -> sparse mask
# =============================================================================

#' Build a GenePathwayMap from a data.frame
#'
#' Converts a two-column (or multi-column) data.frame of gene-pathway
#' associations into a sparse binary mask and a validated \code{GenePathwayMap}
#' S4 object.
#'
#' @param mapping A \code{data.frame} with at least \code{gene_id} and
#'   \code{pathway_id} columns. An optional \code{weight} column (numeric)
#'   will be used as mask values instead of binary 1s.
#' @param feature_genes Optional character vector of gene IDs present in the
#'   expression matrix. If supplied, the mask rows are restricted and ordered
#'   to match these genes exactly. Genes not in the mapping receive all-zero
#'   rows (no pathway membership). This is CRITICAL for aligning the mask
#'   with the columns of a SummarizedExperiment assay.
#' @param min_pathway_size Integer. Pathways with fewer than this many genes
#'   (after intersection with \code{feature_genes}) are dropped. Default 5.
#' @param max_pathway_size Integer. Pathways with more than this many genes
#'   are dropped to reduce noise from overly broad terms. Default 500.
#' @param source Character label for the ontology source. Default "custom".
#'
#' @return A validated \code{GenePathwayMap} object.
#'
#' @details
#' ## The Biology-to-Matrix Pipeline
#'
#' The key idea is straightforward:
#'
#' \preformatted{
#'   Gene-Pathway Table          Sparse Mask (dgCMatrix)
#'   +---------+------------+     Pathways:  P1  P2  P3
#'   | gene_id | pathway_id |     Gene_A  [  1   0   1 ]
#'   +---------+------------+     Gene_B  [  1   1   0 ]
#'   | Gene_A  | P1         |     Gene_C  [  0   1   0 ]
#'   | Gene_A  | P3         |     Gene_D  [  0   0   0 ]  <- not in any pathway
#'   | Gene_B  | P1         |
#'   | Gene_B  | P2         |
#'   | Gene_C  | P2         |
#'   +---------+------------+
#' }
#'
#' Pitfalls this function handles:
#' \enumerate{
#'   \item \strong{ID Mismatch}: Gene IDs in the mapping (e.g., Entrez) may
#'     not match the feature IDs in your expression data (e.g., Ensembl).
#'     Convert beforehand with \code{harmonizeGeneIds()}.
#'   \item \strong{Duplicate Entries}: Many-to-many mappings (one gene in
#'     multiple pathways) are valid and preserved; exact-duplicate rows are
#'     silently collapsed.
#'   \item \strong{Pathway Size Filtering}: Very small or very large pathways
#'     add noise. Filtered by \code{min/max_pathway_size}.
#'   \item \strong{Feature Alignment}: The mask rows MUST correspond 1:1 to
#'     the features (rows) of the expression matrix. The
#'     \code{feature_genes} argument enforces this alignment.
#' }
#'
#' @export
#' @examples
#' mapping_df <- data.frame(
#'     gene_id    = c("ENSG001", "ENSG001", "ENSG002", "ENSG003"),
#'     pathway_id = c("KEGG_P53", "KEGG_APOPTOSIS", "KEGG_P53", "KEGG_APOPTOSIS")
#' )
#' features <- c("ENSG001", "ENSG002", "ENSG003", "ENSG004")
#' gpm <- buildGenePathwayMap(mapping_df, feature_genes = features,
#'     min_pathway_size = 1)
buildGenePathwayMap <- function(mapping,
                                feature_genes = NULL,
                                min_pathway_size = 5L,
                                max_pathway_size = 500L,
                                source = "custom") {

    ## ---- input validation ---------------------------------------------------
    stopifnot(is.data.frame(mapping))
    required <- c("gene_id", "pathway_id")
    missing_cols <- setdiff(required, colnames(mapping))
    if (length(missing_cols) > 0) {
        stop("mapping is missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }

    ## Coerce to character to prevent factor-level mismatches
    mapping$gene_id    <- as.character(mapping$gene_id)
    mapping$pathway_id <- as.character(mapping$pathway_id)

    ## Remove exact-duplicate rows
    mapping <- unique(mapping[, c("gene_id", "pathway_id",
                                  if ("weight" %in% colnames(mapping))
                                      "weight"), drop = FALSE])

    ## ---- intersect with feature space ---------------------------------------
    if (!is.null(feature_genes)) {
        feature_genes <- as.character(feature_genes)
        mapping <- mapping[mapping$gene_id %in% feature_genes, , drop = FALSE]
    }

    if (nrow(mapping) == 0) {
        stop("No gene-pathway pairs remain after intersecting with ",
             "feature_genes. Check that ID types match (Ensembl vs Entrez?).")
    }

    ## ---- pathway size filtering ---------------------------------------------
    pw_sizes <- table(mapping$pathway_id)
    keep_pw <- names(pw_sizes)[pw_sizes >= min_pathway_size &
                                pw_sizes <= max_pathway_size]

    if (length(keep_pw) == 0) {
        stop("No pathways remain after size filtering ",
             "(min=", min_pathway_size, ", max=", max_pathway_size, "). ",
             "Consider relaxing thresholds.")
    }

    mapping <- mapping[mapping$pathway_id %in% keep_pw, , drop = FALSE]
    n_dropped <- length(pw_sizes) - length(keep_pw)
    if (n_dropped > 0) {
        message("Dropped ", n_dropped, " pathway(s) outside size range [",
                min_pathway_size, ", ", max_pathway_size, "].")
    }

    ## ---- build ordered index vectors ----------------------------------------
    if (!is.null(feature_genes)) {
        ## Preserve the exact order from the expression matrix
        all_genes <- feature_genes
    } else {
        all_genes <- sort(unique(mapping$gene_id))
    }
    all_pathways <- sort(unique(mapping$pathway_id))

    gene_idx <- setNames(seq_along(all_genes), all_genes)
    pw_idx   <- setNames(seq_along(all_pathways), all_pathways)

    ## ---- assemble sparse matrix ---------------------------------------------
    row_i <- gene_idx[mapping$gene_id]
    col_j <- pw_idx[mapping$pathway_id]

    ## Use weights if provided, else binary
    vals <- if ("weight" %in% colnames(mapping)) {
        as.numeric(mapping$weight)
    } else {
        rep(1, nrow(mapping))
    }

    mask <- sparseMatrix(
        i    = as.integer(row_i),
        j    = as.integer(col_j),
        x    = vals,
        dims = c(length(all_genes), length(all_pathways)),
        dimnames = list(all_genes, all_pathways)
    )

    ## ---- return validated object ---------------------------------------------
    obj <- new("GenePathwayMap",
        mapping       = mapping,
        mask          = mask,
        gene_index    = setNames(all_genes, names(gene_idx)),
        pathway_index = setNames(all_pathways, names(pw_idx)),
        source        = source
    )
    validObject(obj)
    obj
}


# =============================================================================
# Convenience: MSigDB via msigdbr
# =============================================================================

#' Build a GenePathwayMap from MSigDB
#'
#' Fetches gene sets from the \pkg{msigdbr} package and produces a
#' \code{GenePathwayMap} using Ensembl gene IDs by default.
#'
#' @param species Character. Species name for msigdbr (default "Homo sapiens").
#' @param category Character. MSigDB category (e.g., "H" for Hallmarks,
#'   "C2" for curated, "C5" for GO).
#' @param subcategory Character or NULL. E.g., "CP:KEGG", "GO:BP".
#' @param gene_id_type One of "ensembl_gene", "entrez_gene", "gene_symbol".
#' @param feature_genes Optional character vector for alignment. See
#'   \code{\link{buildGenePathwayMap}}.
#' @inheritParams buildGenePathwayMap
#'
#' @return A \code{GenePathwayMap} object.
#' @export
buildMapFromMSigDB <- function(species = "Homo sapiens",
                                category = "H",
                                subcategory = NULL,
                                gene_id_type = c("ensembl_gene",
                                                 "entrez_gene",
                                                 "gene_symbol"),
                                feature_genes = NULL,
                                min_pathway_size = 5L,
                                max_pathway_size = 500L) {

    gene_id_type <- match.arg(gene_id_type)

    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Package 'msigdbr' is required. Install with:\n",
             "  BiocManager::install('msigdbr')")
    }

    msig_args <- list(species = species, category = category)
    if (!is.null(subcategory)) msig_args$subcategory <- subcategory
    msig_df <- do.call(msigdbr::msigdbr, msig_args)

    ## Map msigdbr column names to our standard
    id_col <- switch(gene_id_type,
        ensembl_gene = "ensembl_gene",
        entrez_gene  = "entrez_gene",
        gene_symbol  = "gene_symbol"
    )

    mapping <- data.frame(
        gene_id    = as.character(msig_df[[id_col]]),
        pathway_id = as.character(msig_df$gs_name),
        stringsAsFactors = FALSE
    )

    src_label <- paste0("MSigDB_", category,
                        if (!is.null(subcategory)) paste0(":", subcategory)
                        else "")

    buildGenePathwayMap(mapping,
        feature_genes    = feature_genes,
        min_pathway_size = min_pathway_size,
        max_pathway_size = max_pathway_size,
        source           = src_label
    )
}


# =============================================================================
# ID Harmonization helper
# =============================================================================

#' Harmonize Gene ID Types
#'
#' Translates between Ensembl, Entrez, and Symbol using
#' \code{org.Hs.eg.db}. Handles the 1:many and missing-ID problems that
#' plague biology-to-matrix pipelines.
#'
#' @param ids Character vector of gene IDs.
#' @param from One of "ENSEMBL", "ENTREZID", "SYMBOL".
#' @param to One of "ENSEMBL", "ENTREZID", "SYMBOL".
#' @param organism OrgDb object. Default \code{org.Hs.eg.db::org.Hs.eg.db}.
#'
#' @return A data.frame with columns \code{from} and \code{to}. Rows with no
#'   mapping are dropped and a message is printed.
#' @export
harmonizeGeneIds <- function(ids,
                              from = "ENSEMBL",
                              to = "SYMBOL",
                              organism = NULL) {

    if (is.null(organism)) {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            stop("Install org.Hs.eg.db:  BiocManager::install('org.Hs.eg.db')")
        }
        organism <- org.Hs.eg.db::org.Hs.eg.db
    }

    result <- AnnotationDbi::select(
        organism,
        keys     = unique(as.character(ids)),
        keytype  = from,
        columns  = to
    )

    ## Report drop rate
    n_input   <- length(unique(ids))
    n_mapped  <- length(unique(result[[from]]))
    n_dropped <- n_input - n_mapped
    if (n_dropped > 0) {
        message(n_dropped, " / ", n_input, " IDs (",
                round(100 * n_dropped / n_input, 1),
                "%) could not be mapped from ", from, " -> ", to, ".")
    }

    colnames(result) <- c("from", "to")
    result[!is.na(result$to), ]
}
