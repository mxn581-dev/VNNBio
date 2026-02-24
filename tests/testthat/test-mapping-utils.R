# ================================================================
# test-mapping-utils.R — Sparse Mask Construction
# These tests require NO Julia installation.
# They cover the most bug-prone part of the whole library.
# ================================================================

library(Matrix)

# ── Helper: reusable test fixtures ───────────────────────────────

make_test_mapping <- function() {
    data.frame(
        gene_id    = c("G1","G1","G2","G2","G3","G4","G5",
                        "G6","G7","G8","G9","G10"),
        pathway_id = c("PA","PB","PA","PB","PA","PA","PA",
                        "PB","PB","PB","PB","PB"),
        stringsAsFactors = FALSE
    )
}

make_large_mapping <- function(n_genes = 100, n_pathways = 10,
                                genes_per_pathway = 15) {
    rows <- lapply(seq_len(n_pathways), function(p) {
        gene_ids <- paste0("Gene_", sample(n_genes, genes_per_pathway))
        data.frame(
            gene_id    = gene_ids,
            pathway_id = paste0("Pathway_", p),
            stringsAsFactors = FALSE
        )
    })
    unique(do.call(rbind, rows))
}


# ── Basic construction ───────────────────────────────────────────

test_that("buildGenePathwayMap creates valid object from simple input", {
    mapping <- make_test_mapping()
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 1)

    expect_s4_class(gpm, "GenePathwayMap")
    expect_true(inherits(gpm@mask, "dgCMatrix"))
    expect_equal(gpm@source, "custom")

    # All genes should be in the mask
    expect_true(all(unique(mapping$gene_id) %in% rownames(gpm@mask)))

    # All pathways should be columns
    expect_true(all(unique(mapping$pathway_id) %in% colnames(gpm@mask)))
})

test_that("mask values are binary (0 or 1) by default", {
    mapping <- make_test_mapping()
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 1)

    # All nonzero values should be exactly 1
    expect_true(all(gpm@mask@x == 1))
})

test_that("mask dimensions match gene and pathway counts", {
    mapping <- make_test_mapping()
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 1)

    n_genes <- length(unique(mapping$gene_id))
    n_pathways <- length(unique(mapping$pathway_id))

    expect_equal(nrow(gpm@mask), n_genes)
    expect_equal(ncol(gpm@mask), n_pathways)
    expect_equal(length(gpm@gene_index), n_genes)
    expect_equal(length(gpm@pathway_index), n_pathways)
})


# ── Feature alignment (THE critical behavior) ────────────────────

test_that("feature_genes forces exact row order and count", {
    mapping <- data.frame(
        gene_id    = c("G1", "G2", "G3", "G1", "G3"),
        pathway_id = c("P1", "P1", "P1", "P2", "P2")
    )
    features <- c("G3", "G1", "G2", "G4", "G5")  # different order + extras

    gpm <- buildGenePathwayMap(mapping,
        feature_genes = features,
        min_pathway_size = 1
    )

    # Mask should have exactly 5 rows (matching features), not 3
    expect_equal(nrow(gpm@mask), 5)

    # Row order must match feature_genes exactly
    expect_equal(rownames(gpm@mask), features)

    # G4 and G5 should have all-zero rows (no pathway membership)
    expect_equal(sum(gpm@mask["G4", ]), 0)
    expect_equal(sum(gpm@mask["G5", ]), 0)

    # G1 should be in both pathways
    expect_equal(sum(gpm@mask["G1", ]), 2)
})

test_that("feature_genes with zero overlap throws informative error", {
    mapping <- data.frame(
        gene_id    = c("ENSG001", "ENSG002"),
        pathway_id = c("P1", "P1")
    )
    features <- c("GENE_A", "GENE_B")  # completely different ID namespace

    expect_error(
        buildGenePathwayMap(mapping, feature_genes = features),
        "No gene-pathway pairs remain"
    )
})

test_that("feature_genes preserves SummarizedExperiment alignment", {
    # Simulate what happens in the real pipeline:
    # rownames(se) -> feature_genes -> mask rows -> Julia input columns
    set.seed(1)
    se_genes <- paste0("Gene_", sprintf("%04d", sample(1:500, 50)))

    mapping <- data.frame(
        gene_id    = rep(se_genes[1:30], each = 2),
        pathway_id = paste0("PW_", rep(1:6, 10))
    )

    gpm <- buildGenePathwayMap(mapping,
        feature_genes = se_genes,
        min_pathway_size = 5
    )

    # THE critical invariant: mask row order == SE row order
    expect_identical(rownames(gpm@mask), se_genes)
    expect_equal(nrow(gpm@mask), length(se_genes))
})


# ── Pathway size filtering ───────────────────────────────────────

test_that("min_pathway_size filters small pathways", {
    mapping <- data.frame(
        gene_id    = c("G1","G2","G3",                 # PW_big has 3 genes
                        "G4"),                          # PW_tiny has 1 gene
        pathway_id = c("PW_big","PW_big","PW_big",
                        "PW_tiny")
    )

    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 2)

    expect_equal(ncol(gpm@mask), 1)  # only PW_big survives
    expect_equal(colnames(gpm@mask), "PW_big")
})

test_that("max_pathway_size filters large pathways", {
    # Pathway with 20 genes
    mapping <- data.frame(
        gene_id    = paste0("G", 1:20),
        pathway_id = rep("PW_huge", 20)
    )

    expect_error(
        buildGenePathwayMap(mapping, min_pathway_size = 1, max_pathway_size = 10),
        "No pathways remain"
    )
})

test_that("size filtering reports dropped count", {
    mapping <- data.frame(
        gene_id    = c(paste0("G", 1:10), "G99"),
        pathway_id = c(rep("PW_ok", 10), "PW_tiny")
    )

    expect_message(
        buildGenePathwayMap(mapping, min_pathway_size = 5),
        "Dropped 1 pathway"
    )
})

test_that("all pathways filtered produces clear error", {
    mapping <- data.frame(
        gene_id    = c("G1", "G2"),
        pathway_id = c("P1", "P1")
    )

    expect_error(
        buildGenePathwayMap(mapping, min_pathway_size = 10),
        "No pathways remain"
    )
})


# ── Deduplication ────────────────────────────────────────────────

test_that("duplicate gene-pathway pairs are collapsed", {
    mapping <- data.frame(
        gene_id    = c("G1","G1","G1","G2","G2","G3","G4","G5"),
        pathway_id = c("P1","P1","P1","P1","P1","P1","P1","P1")  # G1-P1 x3
    )

    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 1)

    # G1 should appear once in the mask, not three times
    expect_equal(gpm@mask["G1", "P1"], 1)
    expect_equal(Matrix::nnzero(gpm@mask), 5)  # 5 unique genes
})


# ── Weighted masks ───────────────────────────────────────────────

test_that("weight column is used when present", {
    mapping <- data.frame(
        gene_id    = c("G1","G2","G3","G4","G5"),
        pathway_id = rep("P1", 5),
        weight     = c(0.5, 1.0, 0.3, 0.8, 0.9)
    )

    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 1)

    expect_equal(gpm@mask["G1", "P1"], 0.5)
    expect_equal(gpm@mask["G2", "P1"], 1.0)
})


# ── Sparse storage efficiency ────────────────────────────────────

test_that("large masks are stored sparsely", {
    mapping <- make_large_mapping(n_genes = 500, n_pathways = 50,
                                   genes_per_pathway = 20)

    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 5)

    # Should be dgCMatrix (compressed sparse column)
    expect_true(inherits(gpm@mask, "dgCMatrix"))

    # Density should be well below 50%
    density <- Matrix::nnzero(gpm@mask) / prod(dim(gpm@mask))
    expect_lt(density, 0.5)
})


# ── Input validation ─────────────────────────────────────────────

test_that("non-data.frame input throws error", {
    expect_error(
        buildGenePathwayMap(list(gene_id = "G1", pathway_id = "P1")),
        "is.data.frame"
    )
})

test_that("missing columns throw informative error", {
    expect_error(
        buildGenePathwayMap(data.frame(gene = "G1", pathway = "P1")),
        "missing required columns"
    )
})

test_that("factor columns are handled (coerced to character)", {
    mapping <- data.frame(
        gene_id    = factor(c("G1","G2","G3","G4","G5")),
        pathway_id = factor(rep("P1", 5))
    )

    # Should not error — factors get coerced to character internally
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 1)
    expect_s4_class(gpm, "GenePathwayMap")
})


# ── Overlapping pathways (many-to-many) ──────────────────────────

test_that("genes in multiple pathways are correctly represented", {
    mapping <- data.frame(
        gene_id    = c("G1","G1","G1",     # G1 in 3 pathways
                        "G2","G3","G4","G5","G6",
                        "G7","G8","G9","G10","G11",
                        "G12","G13","G14","G15","G16"),
        pathway_id = c("PA","PB","PC",
                        rep("PA", 5),
                        rep("PB", 5),
                        rep("PC", 5))
    )

    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 5)

    # G1 should have 1s in PA, PB, and PC
    g1_memberships <- as.numeric(gpm@mask["G1", ])
    expect_equal(sum(g1_memberships), 3)
})


# ── Reproducibility ──────────────────────────────────────────────

test_that("same input produces identical mask", {
    mapping <- make_large_mapping(n_genes = 200, n_pathways = 20,
                                   genes_per_pathway = 15)
    features <- paste0("Gene_", 1:200)

    gpm1 <- buildGenePathwayMap(mapping, feature_genes = features,
                                 min_pathway_size = 5)
    gpm2 <- buildGenePathwayMap(mapping, feature_genes = features,
                                 min_pathway_size = 5)

    expect_identical(gpm1@mask, gpm2@mask)
    expect_identical(gpm1@gene_index, gpm2@gene_index)
    expect_identical(gpm1@pathway_index, gpm2@pathway_index)
})
