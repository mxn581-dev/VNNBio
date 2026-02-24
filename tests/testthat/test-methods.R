# ================================================================
# test-methods.R — S4 Method Tests (show, accessors, coercion)
# No Julia required.
# ================================================================

library(Matrix)
library(SummarizedExperiment)

# ── Helper: build a complete GenePathwayMap ──────────────────────

make_gpm <- function() {
    mapping <- data.frame(
        gene_id    = c(paste0("G", 1:10), paste0("G", 5:15)),
        pathway_id = c(rep("PA", 10), rep("PB", 11))
    )
    buildGenePathwayMap(mapping, min_pathway_size = 5)
}


# ── GenePathwayMap accessors ─────────────────────────────────────

test_that("GenePathwayMap accessors return correct types", {
    gpm <- make_gpm()

    expect_s4_class(gpm, "GenePathwayMap")
    expect_true(is.data.frame(mappingTable(gpm)))
    expect_true(inherits(maskMatrix(gpm), "dgCMatrix"))
    expect_type(geneIndex(gpm), "character")
    expect_type(pathwayIndex(gpm), "character")
    expect_type(maskSource(gpm), "character")
    expect_type(nGenes(gpm), "integer")
    expect_type(nPathways(gpm), "integer")
    expect_type(maskDensity(gpm), "double")
})

test_that("GenePathwayMap dim/length methods work", {
    gpm <- make_gpm()

    expect_length(dim(gpm), 2)
    expect_equal(dim(gpm)[1], nGenes(gpm))
    expect_equal(dim(gpm)[2], nPathways(gpm))
    expect_equal(length(gpm), nrow(mappingTable(gpm)))
})

test_that("maskDensity returns value in [0, 1]", {
    gpm <- make_gpm()
    d <- maskDensity(gpm)
    expect_gte(d, 0)
    expect_lte(d, 1)
})

test_that("GenePathwayMap show method runs without error", {
    gpm <- make_gpm()
    expect_output(show(gpm), "GenePathwayMap object")
    expect_output(show(gpm), "Genes:")
    expect_output(show(gpm), "Pathways:")
    expect_output(show(gpm), "density")
})


# ── VNNArchitecture accessors ────────────────────────────────────

test_that("VNNArchitecture accessors return correct values", {
    gpm <- make_gpm()
    arch <- buildArchitecture(gpm, activation = "relu", n_output = 2L)

    expect_s4_class(arch, "VNNArchitecture")
    expect_type(layerMasks(arch), "list")
    expect_length(layerMasks(arch), 1)
    expect_type(layerNames(arch), "character")
    expect_equal(activationFunction(arch), "relu")
    expect_equal(nOutput(arch), 2L)
    expect_equal(nLayers(arch), 1L)
})

test_that("VNNArchitecture show method runs without error", {
    gpm <- make_gpm()
    arch <- buildArchitecture(gpm)
    expect_output(show(arch), "VNNArchitecture object")
    expect_output(show(arch), "Activation:")
    expect_output(show(arch), "topology")
})

test_that("buildArchitecture creates valid architecture from GenePathwayMap", {
    gpm <- make_gpm()
    arch <- buildArchitecture(gpm)

    expect_s4_class(arch, "VNNArchitecture")
    expect_equal(nrow(layerMasks(arch)[[1]]), nGenes(gpm))
    expect_equal(ncol(layerMasks(arch)[[1]]), nPathways(gpm))
    expect_true(validObject(arch))
})

test_that("validateChain returns TRUE for single-layer", {
    gpm <- make_gpm()
    arch <- buildArchitecture(gpm)
    expect_true(validateChain(arch))
})

test_that("validateChain catches dimension mismatch", {
    m1 <- sparseMatrix(i = 1:5, j = rep(1, 5), x = 1, dims = c(10, 3))
    m2 <- sparseMatrix(i = 1:2, j = rep(1, 2), x = 1, dims = c(5, 2))

    # ncol(m1) = 3 but nrow(m2) = 5 -> validity check catches this at construction
    expect_error(
        new("VNNArchitecture",
            layer_masks = list(a = m1, b = m2),
            layer_names = c("L1", "L2"),
            activation = "tanh", n_output = 1L),
        "mismatch"
    )
})

test_that("VNNArchitecture rejects invalid activation", {
    gpm <- make_gpm()
    expect_error(
        new("VNNArchitecture",
            layer_masks = list(m = maskMatrix(gpm)),
            layer_names = "test",
            activation = "banana",
            n_output = 1L),
        "activation"
    )
})


# ── VNNModel accessors ──────────────────────────────────────────

test_that("VNNModel accessors work on default object", {
    model <- new("VNNModel")

    expect_s4_class(architecture(model), "VNNArchitecture")
    expect_equal(nrow(trainingHistory(model)), 0)
    expect_type(importanceScores(model), "list")
    expect_equal(taskType(model), "classification")
    expect_null(inputSE(model))
})

test_that("VNNModel rejects invalid task type", {
    expect_error(
        new("VNNModel", task = "bogus"),
        "task"
    )
})

test_that("importanceTable returns correct format", {
    model <- new("VNNModel",
        importance_scores = list(
            Pathways = c(PA = 5.0, PB = 3.0, PC = 1.0)
        )
    )

    tbl <- importanceTable(model)
    expect_s3_class(tbl, "data.frame")
    expect_true(all(c("layer", "pathway", "importance", "rank") %in%
                    colnames(tbl)))
    expect_equal(nrow(tbl), 3)
    expect_equal(tbl$pathway[1], "PA")
    expect_equal(tbl$rank, 1:3)
})

test_that("importanceTable handles empty scores", {
    model <- new("VNNModel")
    tbl <- importanceTable(model)
    expect_equal(nrow(tbl), 0)
    expect_true(all(c("layer", "pathway", "importance", "rank") %in%
                    colnames(tbl)))
})

test_that("VNNModel show method runs without error", {
    model <- new("VNNModel",
        task = "classification",
        training_history = data.frame(
            epoch = 1:10,
            train_loss = seq(0.7, 0.3, length.out = 10),
            val_loss = seq(0.75, 0.35, length.out = 10)
        ),
        importance_scores = list(
            Pathways = c(Immune = 4.5, CellCycle = 3.2)
        )
    )
    expect_output(show(model), "VNNModel object")
    expect_output(show(model), "classification")
    expect_output(show(model), "Final loss")
    expect_output(show(model), "Immune")
})
