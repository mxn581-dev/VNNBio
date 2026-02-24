# ================================================================
# test-integration-julia.R — Full Pipeline Integration Tests
#
# These tests REQUIRE a working Julia installation.
# They are automatically skipped if Julia is not available.
# ================================================================

library(Matrix)
library(SummarizedExperiment)

# ── Skip if Julia is not available ───────────────────────────────

skip_if_no_julia <- function() {
    julia_ok <- tryCatch({
        requireNamespace("JuliaConnectoR", quietly = TRUE) &&
            JuliaConnectoR::juliaSetupOk()
    }, error = function(e) FALSE)

    if (!julia_ok) {
        skip("Julia is not available — skipping integration tests.")
    }
}

# ── Helper: make a training-ready dataset ────────────────────────

make_training_data <- function(n_genes = 50, n_samples = 100,
                                n_pathways = 5, seed = 42) {
    set.seed(seed)
    gene_names <- paste0("G_", sprintf("%03d", 1:n_genes))
    expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes,
                   dimnames = list(gene_names, paste0("S", 1:n_samples)))

    se <- SummarizedExperiment(
        assays  = list(logcounts = expr),
        colData = DataFrame(
            label = sample(c("Yes", "No"), n_samples, replace = TRUE)
        )
    )

    mapping <- do.call(rbind, lapply(seq_len(n_pathways), function(p) {
        g <- sample(gene_names, 10 + p)
        data.frame(gene_id = g, pathway_id = paste0("PW_", p),
                   stringsAsFactors = FALSE)
    }))
    mapping <- unique(mapping)

    gpm <- buildGenePathwayMap(mapping,
        feature_genes = rownames(se),
        min_pathway_size = 5
    )

    arch <- new("VNNArchitecture",
        layer_masks = list(genes_to_pw = gpm@mask),
        layer_names = "Pathways",
        activation  = "tanh",
        n_output    = 1L
    )

    list(se = se, gpm = gpm, arch = arch)
}


# ── Tests ────────────────────────────────────────────────────────

test_that("Julia backend initializes", {
    skip_if_no_julia()

    result <- initJuliaBackend(force = TRUE)
    expect_true(result)
})

test_that("sparse mask transfers to Julia without error", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()
    mask <- data$gpm@mask

    # Should not error
    expect_no_error(
        VNNBio:::.transferSparseMask(mask, "test_mask")
    )
})

test_that("full training pipeline runs and returns VNNModel", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()

    model <- trainVNN(
        se           = data$se,
        architecture = data$arch,
        label_col    = "label",
        assay_name   = "logcounts",
        task         = "classification",
        epochs       = 10L,       # few epochs just to test pipeline
        learning_rate = 1e-3,
        batch_size    = 16L,
        val_fraction  = 0.2,
        l1_lambda     = 1e-4,
        seed          = 42L,
        verbose       = FALSE
    )

    expect_s4_class(model, "VNNModel")
    expect_equal(model@task, "classification")
})

test_that("training history has correct shape", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()
    n_epochs <- 15L

    model <- trainVNN(
        se = data$se, architecture = data$arch,
        label_col = "label", epochs = n_epochs,
        verbose = FALSE
    )

    hist <- model@training_history
    expect_equal(nrow(hist), n_epochs)
    expect_true("train_loss" %in% colnames(hist))
    expect_true("val_loss" %in% colnames(hist))

    # Losses should be finite positive numbers
    expect_true(all(is.finite(hist$train_loss)))
    expect_true(all(is.finite(hist$val_loss)))
    expect_true(all(hist$train_loss > 0))
})

test_that("training loss decreases over epochs", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()

    model <- trainVNN(
        se = data$se, architecture = data$arch,
        label_col = "label", epochs = 50L,
        learning_rate = 5e-3,  # higher LR to see clear decrease
        verbose = FALSE
    )

    hist <- model@training_history

    # First epoch loss should be higher than last epoch loss
    # (on random data this isn't guaranteed, so use a loose check)
    first_5  <- mean(hist$train_loss[1:5])
    last_5   <- mean(hist$train_loss[46:50])
    expect_lt(last_5, first_5)
})

test_that("importance scores are extracted with correct names", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()

    model <- trainVNN(
        se = data$se, architecture = data$arch,
        label_col = "label", epochs = 10L,
        verbose = FALSE
    )

    scores <- model@importance_scores
    expect_type(scores, "list")
    expect_length(scores, 1)  # one layer
    expect_true("Pathways" %in% names(scores))

    # Score names should match pathway names
    pw_names <- colnames(data$gpm@mask)
    expect_true(all(names(scores[["Pathways"]]) %in% pw_names))

    # Scores should be non-negative (sum of absolute weights)
    expect_true(all(scores[["Pathways"]] >= 0))
})

test_that("importance scores are sorted descending", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()
    model <- trainVNN(
        se = data$se, architecture = data$arch,
        label_col = "label", epochs = 10L,
        verbose = FALSE
    )

    scores <- model@importance_scores[["Pathways"]]
    # Check that scores are sorted descending
    expect_true(all(diff(scores) <= 0))
})

test_that("plotImportance returns a ggplot object", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()
    model <- trainVNN(
        se = data$se, architecture = data$arch,
        label_col = "label", epochs = 10L,
        verbose = FALSE
    )

    p <- plotImportance(model, top_n = 3)
    expect_s3_class(p, "ggplot")
})

test_that("repeated training with same seed gives same results", {
    skip_if_no_julia()
    initJuliaBackend(force = TRUE)

    data <- make_training_data()

    train_args <- list(
        se = data$se, architecture = data$arch,
        label_col = "label", epochs = 10L,
        seed = 99L, verbose = FALSE
    )

    model1 <- do.call(trainVNN, train_args)
    model2 <- do.call(trainVNN, train_args)

    expect_equal(
        model1@training_history$train_loss,
        model2@training_history$train_loss,
        tolerance = 1e-5
    )
})
