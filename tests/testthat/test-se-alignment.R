# ================================================================
# test-se-alignment.R — SummarizedExperiment ↔ Mask Alignment
# Tests the most common source of bugs: making sure the mask rows
# line up with the expression matrix rows.
# No Julia required.
# ================================================================

library(Matrix)
library(SummarizedExperiment)

# ── Helper: make a minimal SE ────────────────────────────────────

make_test_se <- function(n_genes = 50, n_samples = 20, seed = 42) {
    set.seed(seed)
    gene_names <- paste0("Gene_", sprintf("%03d", 1:n_genes))
    expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes,
                   dimnames = list(gene_names, paste0("S", 1:n_samples)))
    SummarizedExperiment(
        assays  = list(logcounts = expr),
        colData = DataFrame(
            condition = sample(c("A", "B"), n_samples, replace = TRUE)
        )
    )
}

make_test_pathway_df <- function(gene_names, n_pathways = 5,
                                  genes_per_pw = 10) {
    rows <- lapply(seq_len(n_pathways), function(p) {
        g <- sample(gene_names, min(genes_per_pw, length(gene_names)))
        data.frame(gene_id = g, pathway_id = paste0("PW_", p),
                   stringsAsFactors = FALSE)
    })
    unique(do.call(rbind, rows))
}


# ── Alignment tests ──────────────────────────────────────────────

test_that("mask rows match SE rows exactly", {
    se <- make_test_se(50, 20)
    mapping <- make_test_pathway_df(rownames(se))

    gpm <- buildGenePathwayMap(mapping,
        feature_genes    = rownames(se),
        min_pathway_size = 5
    )

    # THE invariant
    expect_identical(rownames(gpm@mask), rownames(se))
    expect_equal(nrow(gpm@mask), nrow(se))
})

test_that("mask × expression matrix multiply has correct dimensions", {
    se <- make_test_se(50, 20)
    mapping <- make_test_pathway_df(rownames(se))

    gpm <- buildGenePathwayMap(mapping,
        feature_genes    = rownames(se),
        min_pathway_size = 5
    )

    # Simulate what Julia does: X [samples × genes] %*% mask [genes × pathways]
    X <- t(as.matrix(assay(se, "logcounts")))  # [20 × 50]
    result <- X %*% gpm@mask                    # [20 × n_pathways]

    expect_equal(nrow(result), ncol(se))           # n_samples
    expect_equal(ncol(result), ncol(gpm@mask))     # n_pathways
})

test_that("genes not in any pathway get zero rows in mask", {
    se <- make_test_se(50, 20)

    # Only map the first 20 genes to pathways
    mapping <- data.frame(
        gene_id    = rep(rownames(se)[1:20], each = 2),
        pathway_id = rep(paste0("PW_", 1:5), 8)
    )

    gpm <- buildGenePathwayMap(mapping,
        feature_genes    = rownames(se),
        min_pathway_size = 5
    )

    # Genes 21-50 should have all-zero rows
    unmapped <- rownames(se)[21:50]
    for (g in unmapped) {
        expect_equal(sum(gpm@mask[g, ]), 0,
            info = paste("Gene", g, "should have zero row"))
    }
})

test_that("transposed expression has correct shape for Julia", {
    se <- make_test_se(100, 30)

    X <- t(as.matrix(assay(se, "logcounts")))

    # Julia expects [n_samples × n_genes]
    expect_equal(nrow(X), 30)   # samples
    expect_equal(ncol(X), 100)  # genes
})

test_that("label extraction handles factors correctly", {
    se <- make_test_se(50, 20)

    y <- colData(se)[["condition"]]
    y_numeric <- as.integer(as.factor(y)) - 1L

    expect_true(all(y_numeric %in% c(0L, 1L)))
    expect_length(y_numeric, ncol(se))
})


# ── Edge case: SE with subset of genes ───────────────────────────

test_that("mask works when SE has fewer genes than mapping", {
    # Mapping knows about 200 genes, but SE only has 50
    full_genes <- paste0("Gene_", sprintf("%03d", 1:200))
    se_genes   <- paste0("Gene_", sprintf("%03d", 1:50))

    mapping <- make_test_pathway_df(full_genes, n_pathways = 10,
                                     genes_per_pw = 30)

    gpm <- buildGenePathwayMap(mapping,
        feature_genes    = se_genes,
        min_pathway_size = 5
    )

    expect_equal(nrow(gpm@mask), 50)
    expect_identical(rownames(gpm@mask), se_genes)
})
