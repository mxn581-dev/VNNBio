# ================================================================
# test-mapping-advanced.R — harmonizeGeneIds & buildMapFromMSigDB
#
# These functions depend on Bioconductor annotation packages
# (org.Hs.eg.db, msigdbr) which are in Suggests. Tests skip
# gracefully when packages are unavailable.
# No Julia required.
# ================================================================

library(Matrix)

# ═══════════════════════════════════════════════════════════════════
#                       harmonizeGeneIds
# ═══════════════════════════════════════════════════════════════════

test_that("harmonizeGeneIds requires org.Hs.eg.db or organism arg", {
    # If org.Hs.eg.db is not installed AND no organism is passed,
    # the function should error with an install instruction
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        expect_error(
            harmonizeGeneIds(c("TP53", "BRCA1"), from = "SYMBOL", to = "ENSEMBL"),
            "org.Hs.eg.db"
        )
    }
})

test_that("harmonizeGeneIds SYMBOL -> ENSEMBL returns correct structure", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    result <- harmonizeGeneIds(
        c("TP53", "BRCA1", "EGFR"),
        from = "SYMBOL",
        to   = "ENSEMBL"
    )

    expect_s3_class(result, "data.frame")
    expect_true(all(c("from", "to") %in% colnames(result)))
    expect_gt(nrow(result), 0)

    # Known mappings — these should always exist
    expect_true("TP53" %in% result$from)
    expect_true("BRCA1" %in% result$from)
})

test_that("harmonizeGeneIds ENSEMBL -> SYMBOL works", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    # ENSG00000141510 = TP53
    result <- harmonizeGeneIds(
        "ENSG00000141510",
        from = "ENSEMBL",
        to   = "SYMBOL"
    )

    expect_gt(nrow(result), 0)
    expect_true("TP53" %in% result$to)
})

test_that("harmonizeGeneIds SYMBOL -> ENTREZID works", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    result <- harmonizeGeneIds(
        c("TP53", "BRCA1"),
        from = "SYMBOL",
        to   = "ENTREZID"
    )

    expect_gt(nrow(result), 0)
    # TP53 Entrez ID is 7157
    expect_true("7157" %in% result$to)
})

test_that("harmonizeGeneIds drops unmapped IDs with message", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    # Mix real and fake IDs
    expect_message(
        result <- harmonizeGeneIds(
            c("TP53", "FAKEGENE123", "NOTREAL456"),
            from = "SYMBOL",
            to   = "ENSEMBL"
        ),
        "could not be mapped"
    )

    # Fake genes should not appear in output
    expect_false("FAKEGENE123" %in% result$from)
    expect_false("NOTREAL456" %in% result$from)
    # Real gene should be there
    expect_true("TP53" %in% result$from)
})

test_that("harmonizeGeneIds handles duplicate input IDs", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    # Passing the same ID multiple times should not cause issues
    result <- harmonizeGeneIds(
        c("TP53", "TP53", "BRCA1"),
        from = "SYMBOL",
        to   = "ENSEMBL"
    )

    expect_gt(nrow(result), 0)
    # unique() is called internally on keys, so we should get
    # mappings for 2 unique genes, not 3 rows of input
})

test_that("harmonizeGeneIds with all-fake IDs returns empty df", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    expect_message(
        result <- harmonizeGeneIds(
            c("ZZZZZ999", "FAKE_GENE_42"),
            from = "SYMBOL",
            to   = "ENSEMBL"
        ),
        "could not be mapped"
    )

    expect_equal(nrow(result), 0)
})

test_that("harmonizeGeneIds accepts custom organism argument", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    # Explicitly passing the organism should work the same
    result <- harmonizeGeneIds(
        "TP53",
        from     = "SYMBOL",
        to       = "ENSEMBL",
        organism = org.Hs.eg.db::org.Hs.eg.db
    )

    expect_gt(nrow(result), 0)
})


# ═══════════════════════════════════════════════════════════════════
#                       buildMapFromMSigDB
# ═══════════════════════════════════════════════════════════════════

test_that("buildMapFromMSigDB requires msigdbr package", {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        expect_error(
            buildMapFromMSigDB(category = "H"),
            "msigdbr"
        )
    }
})

test_that("buildMapFromMSigDB Hallmark returns valid GenePathwayMap", {
    skip_if_not_installed("msigdbr")

    gpm <- buildMapFromMSigDB(
        category     = "H",
        gene_id_type = "gene_symbol"
    )

    expect_s4_class(gpm, "GenePathwayMap")
    expect_true(inherits(maskMatrix(gpm), "dgCMatrix"))

    # Hallmark has 50 gene sets
    expect_equal(nPathways(gpm), 50)

    # Should have thousands of genes
    expect_gt(nGenes(gpm), 1000)

    # Source label should reflect the category
    expect_true(grepl("MSigDB_H", maskSource(gpm)))
})

test_that("buildMapFromMSigDB with feature_genes restricts rows", {
    skip_if_not_installed("msigdbr")

    # First get all genes to pick a realistic subset
    full <- buildMapFromMSigDB(category = "H", gene_id_type = "gene_symbol")
    some_genes <- head(geneIndex(full), 200)

    gpm <- buildMapFromMSigDB(
        category      = "H",
        gene_id_type  = "gene_symbol",
        feature_genes = some_genes
    )

    # Rows should match feature_genes count exactly
    expect_equal(nrow(maskMatrix(gpm)), length(some_genes))
    expect_identical(rownames(maskMatrix(gpm)), some_genes)
})

test_that("buildMapFromMSigDB respects min/max pathway size", {
    skip_if_not_installed("msigdbr")

    # Very restrictive size filter — fewer pathways should survive
    gpm_strict <- buildMapFromMSigDB(
        category         = "H",
        gene_id_type     = "gene_symbol",
        min_pathway_size = 100,
        max_pathway_size = 150
    )

    gpm_loose <- buildMapFromMSigDB(
        category         = "H",
        gene_id_type     = "gene_symbol",
        min_pathway_size = 5,
        max_pathway_size = 500
    )

    expect_lte(nPathways(gpm_strict), nPathways(gpm_loose))
})

test_that("buildMapFromMSigDB gene_id_type ensembl works", {
    skip_if_not_installed("msigdbr")

    gpm <- buildMapFromMSigDB(
        category     = "H",
        gene_id_type = "ensembl_gene"
    )

    expect_s4_class(gpm, "GenePathwayMap")
    # Ensembl IDs start with ENSG
    some_genes <- head(geneIndex(gpm), 10)
    expect_true(all(grepl("^ENSG", some_genes)))
})

test_that("buildMapFromMSigDB gene_id_type entrez works", {
    skip_if_not_installed("msigdbr")

    gpm <- buildMapFromMSigDB(
        category     = "H",
        gene_id_type = "entrez_gene"
    )

    expect_s4_class(gpm, "GenePathwayMap")
    # Entrez IDs are numeric strings
    some_genes <- head(geneIndex(gpm), 10)
    expect_true(all(grepl("^\\d+$", some_genes)))
})

test_that("buildMapFromMSigDB with subcategory works", {
    skip_if_not_installed("msigdbr")

    gpm <- buildMapFromMSigDB(
        category     = "C2",
        subcategory  = "CP:KEGG",
        gene_id_type = "gene_symbol"
    )

    expect_s4_class(gpm, "GenePathwayMap")
    expect_gt(nPathways(gpm), 0)

    # Pathway names should contain KEGG prefix
    pw_names <- pathwayIndex(gpm)
    expect_true(all(grepl("^KEGG_", pw_names)))

    # Source label should include subcategory
    expect_true(grepl("CP:KEGG", maskSource(gpm)))
})

test_that("buildMapFromMSigDB mask values are binary", {
    skip_if_not_installed("msigdbr")

    gpm <- buildMapFromMSigDB(
        category     = "H",
        gene_id_type = "gene_symbol"
    )

    expect_true(all(maskMatrix(gpm)@x == 1))
})

test_that("buildMapFromMSigDB mask density is reasonable", {
    skip_if_not_installed("msigdbr")

    gpm <- buildMapFromMSigDB(
        category     = "H",
        gene_id_type = "gene_symbol"
    )

    d <- maskDensity(gpm)
    # Hallmark with ~4500 genes and 50 pathways: density should be low
    expect_gt(d, 0)
    expect_lt(d, 0.1)  # well under 10%
})
