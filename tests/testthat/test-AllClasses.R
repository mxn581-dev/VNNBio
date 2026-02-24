# ================================================================
# test-AllClasses.R — S4 Class Construction & Validation
# These tests require NO Julia installation.
# ================================================================

library(Matrix)

# ── GenePathwayMap ───────────────────────────────────────────────

test_that("GenePathwayMap constructs with valid inputs", {
    mask <- sparseMatrix(
        i = c(1, 1, 2, 3), j = c(1, 2, 1, 2),
        x = 1, dims = c(3, 2),
        dimnames = list(c("G1", "G2", "G3"), c("P1", "P2"))
    )
    mapping <- data.frame(
        gene_id    = c("G1", "G1", "G2", "G3"),
        pathway_id = c("P1", "P2", "P1", "P2")
    )

    gpm <- new("GenePathwayMap",
        mapping       = mapping,
        mask          = mask,
        gene_index    = c("G1", "G2", "G3"),
        pathway_index = c("P1", "P2"),
        source        = "test"
    )

    expect_s4_class(gpm, "GenePathwayMap")
    expect_equal(nrow(gpm@mask), 3)
    expect_equal(ncol(gpm@mask), 2)
    expect_equal(gpm@source, "test")
})

test_that("GenePathwayMap rejects missing required columns", {
    bad_mapping <- data.frame(gene = "G1", pw = "P1")  # wrong column names

    expect_error(
        new("GenePathwayMap",
            mapping    = bad_mapping,
            mask       = new("dgCMatrix"),
            gene_index = character(0),
            pathway_index = character(0),
            source     = "test"
        ),
        "must contain columns"
    )
})

test_that("GenePathwayMap rejects dimension mismatches", {
    mask <- sparseMatrix(i = 1, j = 1, x = 1, dims = c(3, 2))
    mapping <- data.frame(gene_id = "G1", pathway_id = "P1")

    # gene_index has 2 elements but mask has 3 rows
    expect_error(
        new("GenePathwayMap",
            mapping       = mapping,
            mask          = mask,
            gene_index    = c("G1", "G2"),  # should be 3
            pathway_index = c("P1", "P2"),
            source        = "test"
        ),
        "row count must equal"
    )
})

test_that("GenePathwayMap default prototype is valid", {
    gpm <- new("GenePathwayMap")
    expect_s4_class(gpm, "GenePathwayMap")
    expect_equal(nrow(gpm@mapping), 0)
    expect_equal(gpm@source, "custom")
})


# ── VNNArchitecture ──────────────────────────────────────────────

test_that("VNNArchitecture constructs with defaults", {
    arch <- new("VNNArchitecture")
    expect_s4_class(arch, "VNNArchitecture")
    expect_equal(arch@activation, "tanh")
    expect_equal(arch@n_output, 1L)
    expect_length(arch@layer_masks, 0)
})

test_that("VNNArchitecture stores mask list correctly", {
    m1 <- sparseMatrix(i = 1:3, j = c(1,1,2), x = 1, dims = c(5, 2))
    m2 <- sparseMatrix(i = 1:2, j = c(1,1), x = 1, dims = c(2, 1))

    arch <- new("VNNArchitecture",
        layer_masks = list(genes_to_pathways = m1, pathways_to_system = m2),
        layer_names = c("Pathways", "Systems"),
        activation  = "relu",
        n_output    = 1L
    )

    expect_length(arch@layer_masks, 2)
    expect_equal(names(arch@layer_masks), c("genes_to_pathways", "pathways_to_system"))
    expect_equal(arch@activation, "relu")
})


# ── VNNModel ─────────────────────────────────────────────────────

test_that("VNNModel constructs with defaults", {
    model <- new("VNNModel")
    expect_s4_class(model, "VNNModel")
    expect_null(model@julia_model_ref)
    expect_equal(model@task, "classification")
    expect_equal(nrow(model@training_history), 0)
})
