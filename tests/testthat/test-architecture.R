# ================================================================
# test-architecture.R — VNNArchitecture mask chaining logic
# Tests that multi-layer architectures have compatible dimensions.
# No Julia required.
# ================================================================

library(Matrix)

test_that("single-layer architecture has valid dimensions", {
    mask <- sparseMatrix(
        i = c(1,2,3,1,2), j = c(1,1,1,2,2),
        x = 1, dims = c(100, 10),
        dimnames = list(paste0("G", 1:100), paste0("P", 1:10))
    )

    arch <- new("VNNArchitecture",
        layer_masks = list(input_to_pathway = mask),
        layer_names = "Pathways",
        activation  = "tanh",
        n_output    = 1L
    )

    expect_equal(nrow(arch@layer_masks[[1]]), 100)  # input genes
    expect_equal(ncol(arch@layer_masks[[1]]), 10)   # pathway nodes
})

test_that("multi-layer masks chain correctly (cols_L == rows_L+1)", {
    # Layer 1: 100 genes → 10 pathways
    m1 <- sparseMatrix(
        i = rep(1:100, each = 2),
        j = rep(1:10, 20),
        x = 1, dims = c(100, 10)
    )

    # Layer 2: 10 pathways → 3 systems
    m2 <- sparseMatrix(
        i = c(1:10, 1:10, 1:10),
        j = c(rep(1,10), rep(2,10), rep(3,10)),
        x = 1, dims = c(10, 3)
    )

    # The critical invariant: ncol(m1) == nrow(m2)
    expect_equal(ncol(m1), nrow(m2))

    arch <- new("VNNArchitecture",
        layer_masks = list(genes_to_pw = m1, pw_to_sys = m2),
        layer_names = c("Pathways", "Systems"),
        activation  = "relu",
        n_output    = 1L
    )

    expect_length(arch@layer_masks, 2)
})

test_that("dimension mismatch in multi-layer is detectable", {
    m1 <- sparseMatrix(i = 1, j = 1, x = 1, dims = c(100, 10))
    m2 <- sparseMatrix(i = 1, j = 1, x = 1, dims = c(5, 3))  # 5 != 10

    # ncol(m1) = 10 but nrow(m2) = 5 — this should be caught
    # We test the user-side check (not yet in the class validator,
    # but should be checked before training)
    expect_false(ncol(m1) == nrow(m2))
})

test_that("supported activations are accepted", {
    mask <- sparseMatrix(i = 1:5, j = rep(1, 5), x = 1, dims = c(10, 2))

    for (act in c("tanh", "relu", "sigmoid")) {
        arch <- new("VNNArchitecture",
            layer_masks = list(m = mask),
            layer_names = "test",
            activation  = act,
            n_output    = 1L
        )
        expect_equal(arch@activation, act)
    }
})
