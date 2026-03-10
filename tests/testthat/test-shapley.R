# ================================================================
# test-shapley.R — Shapley Pathway Attribution Tests
#
# All tests construct fake VNN parameters in R, so NO Julia required.
# The key insight: we test with known weight matrices where we can
# predict what the Shapley values should be.
# ================================================================

library(Matrix)

# ═══════════════════════════════════════════════════════════════════
#                      Fixture builders
# ═══════════════════════════════════════════════════════════════════

#' Build a minimal VNN parameter set and model for testing.
#'
#' Architecture: 10 genes -> 3 pathways -> 1 output (classification)
#' Mask: each pathway gets ~3-4 genes, non-overlapping by default.
make_shapley_fixture <- function(seed = 42, overlapping = FALSE) {
    set.seed(seed)

    n_genes <- 10
    n_pathways <- 3
    pw_names <- c("Apoptosis", "CellCycle", "Immune")
    gene_names <- paste0("G", seq_len(n_genes))

    ## Build mask: [genes x pathways]
    mask <- matrix(0, nrow = n_genes, ncol = n_pathways,
                   dimnames = list(gene_names, pw_names))
    mask[1:3, 1] <- 1   # Apoptosis: G1-G3
    mask[4:7, 2] <- 1   # CellCycle: G4-G7
    mask[8:10, 3] <- 1  # Immune: G8-G10

    if (overlapping) {
        mask[3, 2] <- 1  # G3 also in CellCycle
        mask[7, 3] <- 1  # G7 also in Immune
    }

    ## Hidden layer weights (pre-mask application)
    W1 <- matrix(rnorm(n_genes * n_pathways, sd = 0.5),
                 nrow = n_genes, ncol = n_pathways,
                 dimnames = list(gene_names, pw_names))

    b1 <- rnorm(n_pathways, sd = 0.1)
    names(b1) <- pw_names

    ## Output layer: [pathways x 1]
    W_out <- matrix(rnorm(n_pathways), nrow = n_pathways, ncol = 1)
    b_out <- rnorm(1, sd = 0.1)

    ## Construct params list (same structure as .extractVNNParams output)
    params <- list(
        layers = list(
            list(weight = W1, bias = b1, mask = mask),
            list(weight = W_out, bias = b_out, mask = NULL)
        ),
        n_layers   = 2L,
        activation = "tanh",
        task       = "classification"
    )

    ## Build matching GenePathwayMap and VNNModel
    mapping_rows <- lapply(seq_len(n_pathways), function(j) {
        genes_in_pw <- gene_names[mask[, j] == 1]
        data.frame(gene_id = genes_in_pw, pathway_id = pw_names[j],
                   stringsAsFactors = FALSE)
    })
    mapping <- do.call(rbind, mapping_rows)

    gpm <- buildGenePathwayMap(mapping,
        feature_genes    = gene_names,
        min_pathway_size = 1
    )
    arch <- buildArchitecture(gpm)

    ## Fake importance scores (not used by Shapley but needed for valid model)
    scores <- setNames(c(3, 2, 1), pw_names)
    layer_label <- paste0(maskSource(gpm), " pathways")

    model <- new("VNNModel",
        architecture      = arch,
        importance_scores = list("custom pathways" = scores),
        task              = "classification"
    )

    ## Generate some fake expression data [samples x genes]
    n_samples <- 20
    X <- matrix(rnorm(n_samples * n_genes), nrow = n_samples, ncol = n_genes,
                dimnames = list(paste0("S", seq_len(n_samples)), gene_names))

    list(params = params, model = model, X = X, mask = mask,
         pw_names = pw_names, gene_names = gene_names)
}

#' Build a fixture where one pathway is a null player (zero output weight).
make_null_player_fixture <- function(seed = 42) {
    fix <- make_shapley_fixture(seed = seed)

    ## Zero out the output weight for pathway 3 (Immune)
    fix$params$layers[[2]]$weight[3, 1] <- 0

    fix
}

#' Build a fixture where two pathways are symmetric (identical weights).
make_symmetric_fixture <- function(seed = 42) {
    fix <- make_shapley_fixture(seed = seed)

    ## Make pathways 1 and 2 have identical hidden-layer weights
    ## and identical output weights
    fix$params$layers[[1]]$weight[, 2] <- fix$params$layers[[1]]$weight[, 1]
    fix$params$layers[[1]]$bias[2] <- fix$params$layers[[1]]$bias[1]
    fix$params$layers[[2]]$weight[2, 1] <- fix$params$layers[[2]]$weight[1, 1]

    ## Also make their masks identical (same genes)
    fix$params$layers[[1]]$mask[, 2] <- fix$params$layers[[1]]$mask[, 1]

    fix
}


# ═══════════════════════════════════════════════════════════════════
#                 Internal: forward pass & precompute
# ═══════════════════════════════════════════════════════════════════

test_that(".forwardPassR produces scalar output in [0,1] for classification", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    pred <- VNNBio:::.forwardPassR(x, fix$params)

    expect_length(pred, 1)
    expect_true(is.finite(pred))
    expect_gte(pred, 0)
    expect_lte(pred, 1)
})

test_that(".forwardPassR with active_pathways subset gives different result", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    full <- VNNBio:::.forwardPassR(x, fix$params)
    partial <- VNNBio:::.forwardPassR(x, fix$params, active_pathways = c(1L))

    # With only 1 of 3 pathways active, prediction should differ
    expect_false(isTRUE(all.equal(full, partial)))
})

test_that(".forwardPassR with no pathways gives baseline", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    baseline <- VNNBio:::.forwardPassR(x, fix$params,
                                        active_pathways = integer(0))

    expect_length(baseline, 1)
    expect_true(is.finite(baseline))
    expect_gte(baseline, 0)
    expect_lte(baseline, 1)
})

test_that(".precomputePathwayActivations returns correct dimensions", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    precomp <- VNNBio:::.precomputePathwayActivations(x, fix$params)

    expect_length(precomp$h, 3)  # 3 pathways
    expect_true(all(is.finite(precomp$h)))
    expect_equal(nrow(precomp$W_out), 3)
    expect_equal(ncol(precomp$W_out), 1)
})

test_that(".evalCoalition full vs empty differs", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    precomp <- VNNBio:::.precomputePathwayActivations(x, fix$params)

    full_pred <- VNNBio:::.evalCoalition(precomp, 1:3)
    empty_pred <- VNNBio:::.evalCoalition(precomp, integer(0))

    expect_false(isTRUE(all.equal(full_pred, empty_pred)))
})

test_that(".evalCoalition matches .forwardPassR", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    ## Full forward pass
    fwd_full <- VNNBio:::.forwardPassR(x, fix$params)

    ## Coalition-based full evaluation
    precomp <- VNNBio:::.precomputePathwayActivations(x, fix$params)
    coal_full <- VNNBio:::.evalCoalition(precomp, 1:3)

    expect_equal(fwd_full, coal_full, tolerance = 1e-10)
})

test_that(".evalCoalition single pathway matches .forwardPassR with that pathway", {
    fix <- make_shapley_fixture()
    x <- fix$X[1, ]

    for (j in 1:3) {
        fwd <- VNNBio:::.forwardPassR(x, fix$params, active_pathways = j)
        precomp <- VNNBio:::.precomputePathwayActivations(x, fix$params)
        coal <- VNNBio:::.evalCoalition(precomp, j)
        expect_equal(fwd, coal, tolerance = 1e-10,
                     info = paste("Pathway", j))
    }
})


# ═══════════════════════════════════════════════════════════════════
#                 shapleyPathwayAttribution
# ═══════════════════════════════════════════════════════════════════

test_that("shapleyPathwayAttribution returns correct structure (global)", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "global",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    expect_type(shap, "list")
    expect_true("global" %in% names(shap))
    expect_true("predictions" %in% names(shap))
    expect_true("baseline" %in% names(shap))
    expect_true("pathway_names" %in% names(shap))

    expect_length(shap$global, 3)
    expect_length(shap$predictions, nrow(fix$X))
    expect_null(shap$local)  # not requested
})

test_that("shapleyPathwayAttribution returns correct structure (local)", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    expect_true("local" %in% names(shap))
    expect_null(shap$global)

    expect_equal(nrow(shap$local), nrow(fix$X))
    expect_equal(ncol(shap$local), 3)
    expect_equal(colnames(shap$local), fix$pw_names)
})

test_that("shapleyPathwayAttribution mode='both' returns global and local", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "both",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    expect_true("global" %in% names(shap))
    expect_true("local" %in% names(shap))
})

test_that("global Shapley values are non-negative (mean |SHAP|)", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "global",
        n_perm = 100, params = fix$params, verbose = FALSE
    )

    expect_true(all(shap$global >= 0))
})

test_that("global Shapley values are sorted descending", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "global",
        n_perm = 100, params = fix$params, verbose = FALSE
    )

    expect_true(all(diff(shap$global) <= 0))
})


# ═══════════════════════════════════════════════════════════════════
#                   Shapley axiom verification
# ═══════════════════════════════════════════════════════════════════

test_that("EFFICIENCY: Shapley values sum to prediction - baseline", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 500, params = fix$params, seed = 42L, verbose = FALSE
    )

    ## For each sample, sum of SHAP values should equal pred - baseline
    shap_sums <- rowSums(shap$local)
    expected_diffs <- shap$predictions - shap$baseline

    ## With 500 permutations this should be quite close
    expect_equal(shap_sums, expected_diffs, tolerance = 0.05)
})

test_that("NULL PLAYER: pathway with zero output weight gets ~zero SHAP", {
    fix <- make_null_player_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 200, params = fix$params, verbose = FALSE
    )

    ## Pathway 3 (Immune) has zero output weight
    immune_shap <- shap$local[, "Immune"]

    ## All Shapley values for Immune should be ~0
    expect_true(all(abs(immune_shap) < 0.01),
                info = paste("Max |SHAP| for null player:",
                             max(abs(immune_shap))))
})

test_that("SYMMETRY: identical pathways get similar SHAP values", {
    fix <- make_symmetric_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 500, params = fix$params, verbose = FALSE
    )

    ## Pathways 1 (Apoptosis) and 2 (CellCycle) have identical params
    pw1_shap <- shap$local[, "Apoptosis"]
    pw2_shap <- shap$local[, "CellCycle"]

    ## They should have very similar Shapley values per sample
    expect_equal(pw1_shap, pw2_shap, tolerance = 0.05)
})


# ═══════════════════════════════════════════════════════════════════
#                    Reproducibility & stability
# ═══════════════════════════════════════════════════════════════════

test_that("same seed gives identical results", {
    fix <- make_shapley_fixture()

    s1 <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 50, params = fix$params, seed = 99L, verbose = FALSE
    )
    s2 <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 50, params = fix$params, seed = 99L, verbose = FALSE
    )

    expect_identical(s1$local, s2$local)
})

test_that("different seeds give different results", {
    fix <- make_shapley_fixture()

    s1 <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 50, params = fix$params, seed = 1L, verbose = FALSE
    )
    s2 <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "local",
        n_perm = 50, params = fix$params, seed = 999L, verbose = FALSE
    )

    expect_false(identical(s1$local, s2$local))
})

test_that("more permutations reduce variance", {
    fix <- make_shapley_fixture()
    x_single <- fix$X[1, , drop = FALSE]

    ## Run multiple times with few perms
    low_results <- replicate(10, {
        shapleyPathwayAttribution(
            fix$model, x_single, mode = "local",
            n_perm = 20, params = fix$params,
            seed = sample(10000, 1), verbose = FALSE
        )$local[1, ]
    })

    ## Run multiple times with many perms
    high_results <- replicate(10, {
        shapleyPathwayAttribution(
            fix$model, x_single, mode = "local",
            n_perm = 500, params = fix$params,
            seed = sample(10000, 1), verbose = FALSE
        )$local[1, ]
    })

    ## Variance across runs should be smaller with more perms
    var_low <- mean(apply(low_results, 1, var))
    var_high <- mean(apply(high_results, 1, var))

    expect_lt(var_high, var_low)
})


# ═══════════════════════════════════════════════════════════════════
#                       Error handling
# ═══════════════════════════════════════════════════════════════════

test_that("n_perm < 10 throws error", {
    fix <- make_shapley_fixture()

    expect_error(
        shapleyPathwayAttribution(
            fix$model, fix$X, n_perm = 5,
            params = fix$params, verbose = FALSE
        ),
        "n_perm must be >= 10"
    )
})

test_that("model without Julia ref or params throws error", {
    fix <- make_shapley_fixture()

    expect_error(
        shapleyPathwayAttribution(
            fix$model, fix$X, verbose = FALSE
        ),
        "no Julia reference"
    )
})

test_that("single-layer params (no output layer) throws error", {
    fix <- make_shapley_fixture()
    fix$params$layers <- fix$params$layers[1]  # remove output layer
    fix$params$n_layers <- 1L

    expect_error(
        shapleyPathwayAttribution(
            fix$model, fix$X, params = fix$params, verbose = FALSE
        ),
        "at least 2 layers"
    )
})


# ═══════════════════════════════════════════════════════════════════
#                   Visualization functions
# ═══════════════════════════════════════════════════════════════════

test_that("plotShapleyGlobal returns ggplot", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "both",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    p <- plotShapleyGlobal(shap)
    expect_s3_class(p, "gg")
})

test_that("plotShapleyGlobal with top_n works", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "both",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    p <- plotShapleyGlobal(shap, top_n = 2)
    expect_s3_class(p, "gg")
})

test_that("plotShapleyGlobal errors without global values", {
    expect_error(
        plotShapleyGlobal(list(local = matrix(1))),
        "No global Shapley"
    )
})

test_that("plotShapleyLocal returns ggplot", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "both",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    p <- plotShapleyLocal(shap, sample_idx = 1)
    expect_s3_class(p, "gg")
})

test_that("plotShapleyLocal with custom label works", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "both",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    p <- plotShapleyLocal(shap, sample_idx = 3, sample_label = "Patient_003")
    expect_s3_class(p, "gg")
})

test_that("plotShapleyLocal errors without local values", {
    expect_error(
        plotShapleyLocal(list(global = c(a = 1))),
        "No local Shapley"
    )
})

test_that("plotShapleySummary returns ggplot", {
    fix <- make_shapley_fixture()

    shap <- shapleyPathwayAttribution(
        fix$model, fix$X, mode = "both",
        n_perm = 50, params = fix$params, verbose = FALSE
    )

    p <- plotShapleySummary(shap)
    expect_s3_class(p, "gg")
})

test_that("plotShapleySummary errors without local values", {
    expect_error(
        plotShapleySummary(list(global = c(a = 1))),
        "No local Shapley"
    )
})
