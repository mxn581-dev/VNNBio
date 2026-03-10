# ================================================================
# test-significance.R — Permutation Significance Testing
#
# All tests use pre-constructed weight matrices (weight_matrix param)
# so NO Julia is required.
# ================================================================

library(Matrix)

# ── Shared fixture builder ───────────────────────────────────────

make_perm_fixture <- function(n_genes = 30, n_pathways = 3,
                               genes_per_pw = 10, seed = 42) {
    set.seed(seed)

    gene_names <- paste0("G", seq_len(n_genes))
    pw_names <- paste0("PW_", LETTERS[seq_len(n_pathways)])

    # Build mapping with overlapping gene sets
    rows <- lapply(seq_len(n_pathways), function(p) {
        start <- (p - 1) * (genes_per_pw - 3) + 1
        end <- min(start + genes_per_pw - 1, n_genes)
        data.frame(
            gene_id    = gene_names[start:end],
            pathway_id = pw_names[p],
            stringsAsFactors = FALSE
        )
    })
    mapping <- unique(do.call(rbind, rows))

    gpm <- buildGenePathwayMap(mapping,
        feature_genes    = gene_names,
        min_pathway_size = 3
    )
    arch <- buildArchitecture(gpm)
    mask <- as.matrix(maskMatrix(gpm))

    # Build a fake weight matrix with one pathway having much higher weights
    # PW_A gets large weights, others get small — simulates real signal
    W <- matrix(0, nrow = nrow(mask), ncol = ncol(mask),
                dimnames = dimnames(mask))
    for (j in seq_len(ncol(mask))) {
        nonzero_rows <- which(mask[, j] != 0)
        if (j == 1) {
            # First pathway: strong signal
            W[nonzero_rows, j] <- runif(length(nonzero_rows), 0.5, 2.0)
        } else {
            # Other pathways: weak weights
            W[nonzero_rows, j] <- runif(length(nonzero_rows), 0.01, 0.1)
        }
    }

    scores <- colSums(abs(W))
    layer_label <- paste0(maskSource(gpm), " pathways")

    model <- new("VNNModel",
        architecture      = arch,
        importance_scores = setNames(list(scores), layer_label),
        task              = "classification"
    )

    list(model = model, W = W, mask = mask, gpm = gpm)
}


# ═══════════════════════════════════════════════════════════════════
#                   permutePathwayImportance
# ═══════════════════════════════════════════════════════════════════

test_that("permutePathwayImportance returns correct data.frame structure", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 100,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    expect_s3_class(result, "data.frame")
    expected_cols <- c("pathway", "observed", "pathway_size", "null_mean",
                       "null_sd", "z_score", "p_value", "p_adj")
    expect_true(all(expected_cols %in% colnames(result)))
})

test_that("permutePathwayImportance returns one row per pathway", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 50,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    n_pw <- ncol(fix$mask)
    expect_equal(nrow(result), n_pw)
})

test_that("permutePathwayImportance results sorted by p_value ascending", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 200,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    expect_true(all(diff(result$p_value) >= 0))
})

test_that("observed importance matches colSums(abs(W))", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 50,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    expected_obs <- colSums(abs(fix$W))
    # Match by pathway name since result is sorted by p_value
    for (pw in names(expected_obs)) {
        row <- result[result$pathway == pw, ]
        expect_equal(row$observed, unname(expected_obs[pw]),
                     tolerance = 1e-10,
                     info = paste("Pathway:", pw))
    }
})

test_that("pathway_size matches mask column sums", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 50,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    expected_sizes <- colSums(fix$mask != 0)
    for (pw in names(expected_sizes)) {
        row <- result[result$pathway == pw, ]
        expect_equal(row$pathway_size, unname(expected_sizes[pw]),
                     info = paste("Pathway:", pw))
    }
})

test_that("p_values are in (0, 1]", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 100,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    # Pseudocount ensures p > 0
    expect_true(all(result$p_value > 0))
    expect_true(all(result$p_value <= 1))
    expect_true(all(result$p_adj > 0))
    expect_true(all(result$p_adj <= 1))
})

test_that("p_adj >= p_value (BH correction only inflates)", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 200,
        weight_matrix = fix$W,
        verbose       = FALSE
    )

    expect_true(all(result$p_adj >= result$p_value - 1e-12))
})

test_that("strong-signal pathway has lower p-value than weak pathways", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model,
        n_perm        = 500,
        weight_matrix = fix$W,
        seed          = 42L,
        verbose       = FALSE
    )

    # PW_A was given large weights (0.5-2.0), others got small (0.01-0.1)
    pw_a <- result[result$pathway == "PW_A", ]
    others <- result[result$pathway != "PW_A", ]

    expect_lt(pw_a$p_value, max(others$p_value))
    expect_gt(pw_a$z_score, min(others$z_score))
})

test_that("reproducible results with same seed", {
    fix <- make_perm_fixture()

    r1 <- permutePathwayImportance(
        fix$model, n_perm = 100, seed = 99L,
        weight_matrix = fix$W, verbose = FALSE
    )
    r2 <- permutePathwayImportance(
        fix$model, n_perm = 100, seed = 99L,
        weight_matrix = fix$W, verbose = FALSE
    )

    expect_identical(r1$p_value, r2$p_value)
    expect_identical(r1$z_score, r2$z_score)
})

test_that("different seeds give different null distributions", {
    fix <- make_perm_fixture()

    r1 <- permutePathwayImportance(
        fix$model, n_perm = 100, seed = 1L,
        weight_matrix = fix$W, verbose = FALSE
    )
    r2 <- permutePathwayImportance(
        fix$model, n_perm = 100, seed = 999L,
        weight_matrix = fix$W, verbose = FALSE
    )

    # Null means should differ slightly (different shuffles)
    expect_false(identical(r1$null_mean, r2$null_mean))
})

test_that("more permutations give more precise p-values", {
    fix <- make_perm_fixture()

    r_small <- permutePathwayImportance(
        fix$model, n_perm = 50, seed = 42L,
        weight_matrix = fix$W, verbose = FALSE
    )
    r_large <- permutePathwayImportance(
        fix$model, n_perm = 5000, seed = 42L,
        weight_matrix = fix$W, verbose = FALSE
    )

    # With more permutations, the minimum achievable p-value is smaller
    # p_min = 1 / (n_perm + 1)
    expect_lt(min(r_large$p_value), min(r_small$p_value) + 0.01)
})

test_that("z_scores are finite", {
    fix <- make_perm_fixture()

    result <- permutePathwayImportance(
        fix$model, n_perm = 100,
        weight_matrix = fix$W, verbose = FALSE
    )

    expect_true(all(is.finite(result$z_score)))
})


# ── Error handling ───────────────────────────────────────────────

test_that("n_perm < 10 throws error", {
    fix <- make_perm_fixture()

    expect_error(
        permutePathwayImportance(fix$model, n_perm = 5,
                                  weight_matrix = fix$W, verbose = FALSE),
        "n_perm must be >= 10"
    )
})

test_that("invalid layer index throws error", {
    fix <- make_perm_fixture()

    expect_error(
        permutePathwayImportance(fix$model, layer = 99L,
                                  weight_matrix = fix$W, verbose = FALSE),
        "only has"
    )
})

test_that("dimension mismatch between W and mask throws error", {
    fix <- make_perm_fixture()
    W_wrong <- matrix(0, nrow = 5, ncol = 2)

    expect_error(
        permutePathwayImportance(fix$model, weight_matrix = W_wrong,
                                  verbose = FALSE),
        "dimensions.*do not match"
    )
})

test_that("model without Julia ref or weight_matrix throws error", {
    fix <- make_perm_fixture()

    expect_error(
        permutePathwayImportance(fix$model, verbose = FALSE),
        "no Julia reference"
    )
})


# ═══════════════════════════════════════════════════════════════════
#                    plotPermutationResults
# ═══════════════════════════════════════════════════════════════════

test_that("plotPermutationResults returns ggplot", {
    fix <- make_perm_fixture()

    perm <- permutePathwayImportance(
        fix$model, n_perm = 100,
        weight_matrix = fix$W, verbose = FALSE
    )

    p <- plotPermutationResults(perm)
    expect_s3_class(p, "gg")
})

test_that("plotPermutationResults respects top_n", {
    fix <- make_perm_fixture()

    perm <- permutePathwayImportance(
        fix$model, n_perm = 100,
        weight_matrix = fix$W, verbose = FALSE
    )

    p <- plotPermutationResults(perm, top_n = 1)
    expect_s3_class(p, "gg")

    # Extract plot data — should only have 1 pathway
    plot_data <- ggplot2::ggplot_build(p)$data[[1]]
    expect_equal(nrow(plot_data), 1)
})

test_that("plotPermutationResults with show_null=FALSE works", {
    fix <- make_perm_fixture()

    perm <- permutePathwayImportance(
        fix$model, n_perm = 100,
        weight_matrix = fix$W, verbose = FALSE
    )

    p <- plotPermutationResults(perm, show_null = FALSE)
    expect_s3_class(p, "gg")
})

test_that("plotPermutationResults errors on bad input", {
    expect_error(
        plotPermutationResults(data.frame(x = 1)),
        "missing required columns"
    )
})

test_that("plotPermutationResults custom alpha works", {
    fix <- make_perm_fixture()

    perm <- permutePathwayImportance(
        fix$model, n_perm = 100,
        weight_matrix = fix$W, verbose = FALSE
    )

    p <- plotPermutationResults(perm, alpha = 0.01)
    expect_s3_class(p, "gg")
})
