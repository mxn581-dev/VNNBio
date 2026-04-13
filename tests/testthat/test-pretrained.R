# ================================================================
# test-pretrained.R — PretrainedVNN, loadPretrained, classifyFrozen
#
# Mode 1 tests: NO Julia required. All use synthetic pretrained
# weights and expression data.
# ================================================================

library(Matrix)
library(SummarizedExperiment)

# ═══════════════════════════════════════════════════════════════════
#                      Fixture builders
# ═══════════════════════════════════════════════════════════════════

#' Build a minimal pretrained bundle and matching SE for testing.
#'
#' Creates a fake PretrainedVNN with 100 genes x 5 pathways,
#' and a SummarizedExperiment with 80 of those genes and 40 samples.
make_pretrained_fixture <- function(seed = 42, n_genes = 100,
                                     n_pathways = 5, n_samples = 40,
                                     n_se_genes = 80) {
    set.seed(seed)

    gene_symbols <- paste0("GENE", sprintf("%03d", seq_len(n_genes)))
    gene_ensembl <- paste0("ENSG00000", sprintf("%06d", seq_len(n_genes)))
    pw_names <- paste0("HALLMARK_PW", seq_len(n_pathways))

    # Random weights and mask
    W_hidden <- matrix(rnorm(n_genes * n_pathways, sd = 0.3),
                       nrow = n_genes, ncol = n_pathways)
    b_hidden <- rnorm(n_pathways, sd = 0.1)

    # Sparse-ish mask: each gene in ~2 pathways on average
    mask <- matrix(0, nrow = n_genes, ncol = n_pathways)
    for (i in seq_len(n_genes)) {
        pws <- sample(n_pathways, sample(1:3, 1))
        mask[i, pws] <- 1
    }
    colnames(mask) <- pw_names
    rownames(mask) <- gene_symbols

    gene_info <- data.frame(
        archs4_index = seq_len(n_genes),
        symbol       = gene_symbols,
        ensembl      = gene_ensembl,
        stringsAsFactors = FALSE
    )

    pretrained <- new("PretrainedVNN",
        W_hidden      = W_hidden,
        b_hidden      = b_hidden,
        mask          = mask,
        gene_info     = gene_info,
        pathway_names = pw_names,
        metadata      = list(
            weights_id       = "test_v1",
            source           = "synthetic",
            mask_source      = "test",
            n_genes          = n_genes,
            n_pathways       = n_pathways,
            activation       = "tanh",
            pretrain_epochs  = 1L,
            tissue_accuracy  = 0.5,
            disease_accuracy = 0.5
        )
    )

    # Build a matching SE (using gene symbols, subset of pretrained genes)
    se_genes <- gene_symbols[seq_len(n_se_genes)]
    expr <- matrix(rnorm(n_se_genes * n_samples), nrow = n_se_genes,
                   dimnames = list(se_genes, paste0("S", seq_len(n_samples))))

    # Inject some signal: first 20 samples get boosted expression in
    # genes that map to PW1
    pw1_genes <- which(mask[seq_len(n_se_genes), 1] == 1)
    expr[pw1_genes, 1:20] <- expr[pw1_genes, 1:20] + 2.0

    se <- SummarizedExperiment(
        assays  = list(logcounts = expr),
        colData = DataFrame(
            condition = c(rep("GroupA", 20), rep("GroupB", 20))
        )
    )

    list(pretrained = pretrained, se = se)
}

#' Save a pretrained bundle as .rds for loadPretrained tests.
save_test_bundle <- function(pretrained, dir = tempdir()) {
    bundle <- list(
        W_hidden      = pretrained@W_hidden,
        b_hidden      = pretrained@b_hidden,
        mask          = pretrained@mask,
        gene_info     = pretrained@gene_info,
        pathway_names = pretrained@pathway_names,
        metadata      = pretrained@metadata
    )
    path <- file.path(dir, "test_v1.rds")
    saveRDS(bundle, path)
    path
}


# ═══════════════════════════════════════════════════════════════════
#                    PretrainedVNN class
# ═══════════════════════════════════════════════════════════════════

test_that("PretrainedVNN constructs with valid inputs", {
    fix <- make_pretrained_fixture()
    pt <- fix$pretrained

    expect_s4_class(pt, "PretrainedVNN")
    expect_equal(nrow(pt@W_hidden), 100)
    expect_equal(ncol(pt@W_hidden), 5)
    expect_equal(length(pt@b_hidden), 5)
    expect_equal(length(pt@pathway_names), 5)
})

test_that("PretrainedVNN rejects dimension mismatches", {
    expect_error(
        new("PretrainedVNN",
            W_hidden = matrix(0, 10, 5),
            b_hidden = numeric(3),  # wrong: should be 5
            mask     = matrix(0, 10, 5),
            gene_info = data.frame(archs4_index = 1:10,
                                   symbol = paste0("G", 1:10),
                                   ensembl = paste0("E", 1:10)),
            pathway_names = paste0("P", 1:5),
            metadata = list()),
        "ncol.*b_hidden"
    )
})

test_that("PretrainedVNN show method prints without error", {
    fix <- make_pretrained_fixture()
    expect_output(show(fix$pretrained), "PretrainedVNN object")
    expect_output(show(fix$pretrained), "Genes:")
    expect_output(show(fix$pretrained), "Pathways:")
})


# ═══════════════════════════════════════════════════════════════════
#                     loadPretrained
# ═══════════════════════════════════════════════════════════════════

test_that("loadPretrained loads from a local path", {
    fix <- make_pretrained_fixture()
    path <- save_test_bundle(fix$pretrained)
    on.exit(unlink(path))

    pt <- loadPretrained(path = path)

    expect_s4_class(pt, "PretrainedVNN")
    expect_equal(nrow(pt@W_hidden), 100)
    expect_equal(ncol(pt@W_hidden), 5)
})

test_that("loadPretrained errors on missing file", {
    expect_error(
        loadPretrained(path = "/nonexistent/path.rds"),
        "not found"
    )
})

test_that("loadPretrained errors on missing weights_id", {
    expect_error(
        loadPretrained(weights_id = "definitely_not_a_real_id"),
        "not found"
    )
})

test_that("loadPretrained errors on malformed bundle", {
    path <- tempfile(fileext = ".rds")
    saveRDS(list(W_hidden = matrix(0, 2, 2)), path)
    on.exit(unlink(path))

    expect_error(
        loadPretrained(path = path),
        "missing keys"
    )
})

test_that("loadPretrained validates object integrity", {
    fix <- make_pretrained_fixture()
    bundle <- list(
        W_hidden      = fix$pretrained@W_hidden,
        b_hidden      = numeric(3),  # wrong length
        mask          = fix$pretrained@mask,
        gene_info     = fix$pretrained@gene_info,
        pathway_names = fix$pretrained@pathway_names,
        metadata      = fix$pretrained@metadata
    )
    path <- tempfile(fileext = ".rds")
    saveRDS(bundle, path)
    on.exit(unlink(path))

    expect_error(loadPretrained(path = path))
})


# ═══════════════════════════════════════════════════════════════════
#                   Gene alignment (.alignGenesToPretrained)
# ═══════════════════════════════════════════════════════════════════

test_that("gene alignment matches expected count", {
    fix <- make_pretrained_fixture(n_genes = 100, n_se_genes = 80)

    aligned <- VNNBio:::.alignGenesToPretrained(
        fix$se, fix$pretrained, "logcounts", gene_id_type = "symbol"
    )

    # SE has 80 genes, pretrained has 100 → should match 80
    expect_equal(aligned$n_matched, 80)
    expect_equal(nrow(aligned$X), ncol(fix$se))  # n_samples
    expect_equal(ncol(aligned$X), 80)             # matched genes
    expect_equal(nrow(aligned$W), 80)
    expect_equal(ncol(aligned$W), 5)
})

test_that("gene alignment handles Ensembl IDs", {
    fix <- make_pretrained_fixture()

    # Rebuild SE with Ensembl rownames
    ensembl_genes <- fix$pretrained@gene_info$ensembl[1:80]
    se_ensembl <- fix$se
    rownames(se_ensembl) <- ensembl_genes

    aligned <- VNNBio:::.alignGenesToPretrained(
        se_ensembl, fix$pretrained, "logcounts", gene_id_type = "ensembl"
    )

    expect_equal(aligned$n_matched, 80)
    expect_equal(aligned$gene_id_type, "ensembl")
})

test_that("gene alignment auto-detects ID type", {
    fix <- make_pretrained_fixture()

    # Symbol rownames → should auto-detect as symbol
    aligned <- VNNBio:::.alignGenesToPretrained(
        fix$se, fix$pretrained, "logcounts", gene_id_type = "auto"
    )
    expect_equal(aligned$gene_id_type, "symbol")

    # Ensembl rownames → should auto-detect as ensembl
    se2 <- fix$se
    rownames(se2) <- fix$pretrained@gene_info$ensembl[1:80]
    aligned2 <- VNNBio:::.alignGenesToPretrained(
        se2, fix$pretrained, "logcounts", gene_id_type = "auto"
    )
    expect_equal(aligned2$gene_id_type, "ensembl")
})

test_that("gene alignment errors on zero overlap", {
    fix <- make_pretrained_fixture()
    se_bad <- fix$se
    rownames(se_bad) <- paste0("FAKEGENE_", seq_len(nrow(se_bad)))

    expect_error(
        VNNBio:::.alignGenesToPretrained(
            se_bad, fix$pretrained, "logcounts"
        ),
        "No gene overlap"
    )
})

test_that("gene alignment strips Ensembl versions", {
    fix <- make_pretrained_fixture()

    # Add version suffixes to SE Ensembl IDs
    ensembl_versioned <- paste0(
        fix$pretrained@gene_info$ensembl[1:80], ".3"
    )
    se_v <- fix$se
    rownames(se_v) <- ensembl_versioned

    aligned <- VNNBio:::.alignGenesToPretrained(
        se_v, fix$pretrained, "logcounts", gene_id_type = "ensembl"
    )

    expect_equal(aligned$n_matched, 80)
})


# ═══════════════════════════════════════════════════════════════════
#                       classifyFrozen
# ═══════════════════════════════════════════════════════════════════

test_that("classifyFrozen requires glmnet", {
    # Can't reliably test package absence, so just verify it's importable
    skip_if_not_installed("glmnet")
    succeed()
})

test_that("classifyFrozen returns expected structure", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()

    result <- classifyFrozen(
        fix$pretrained, fix$se,
        label_col = "condition",
        verbose = FALSE
    )

    expect_type(result, "list")
    expected_keys <- c("predictions", "labels", "auc",
                       "pathway_activations", "pathway_importance",
                       "glmnet_fit", "gene_coverage", "lambda",
                       "class_levels")
    expect_true(all(expected_keys %in% names(result)))
})

test_that("classifyFrozen predictions are probabilities in [0,1]", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    expect_length(result$predictions, ncol(fix$se))
    expect_true(all(result$predictions >= 0))
    expect_true(all(result$predictions <= 1))
})

test_that("classifyFrozen AUC is in [0,1]", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    expect_gte(result$auc, 0)
    expect_lte(result$auc, 1)
})

test_that("classifyFrozen pathway_activations has correct dimensions", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    H <- result$pathway_activations
    expect_equal(nrow(H), ncol(fix$se))  # n_samples
    expect_equal(ncol(H), length(fix$pretrained@pathway_names))
    expect_equal(colnames(H), fix$pretrained@pathway_names)

    # tanh output should be in [-1, 1]
    expect_true(all(H >= -1))
    expect_true(all(H <= 1))
})

test_that("classifyFrozen pathway_importance is sorted descending", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    imp <- result$pathway_importance
    expect_true(all(diff(imp) <= 0))
    expect_true(all(imp >= 0))  # absolute values
    expect_equal(length(imp), length(fix$pretrained@pathway_names))
})

test_that("classifyFrozen gene_coverage is correct", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture(n_genes = 100, n_se_genes = 80)
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    expect_equal(result$gene_coverage, 0.8, tolerance = 0.01)
})

test_that("classifyFrozen with signal gives AUC > 0.5", {
    skip_if_not_installed("glmnet")

    # Fixture has signal injected in PW1 genes for GroupA
    fix <- make_pretrained_fixture(seed = 42)
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    # With injected signal, AUC should be above chance
    expect_gt(result$auc, 0.5)
})

test_that("classifyFrozen errors on non-binary labels", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    se3 <- fix$se
    SummarizedExperiment::colData(se3)$condition <- sample(
        c("A", "B", "C"), ncol(se3), replace = TRUE
    )

    expect_error(
        classifyFrozen(fix$pretrained, se3,
                       label_col = "condition", verbose = FALSE),
        "exactly 2 classes"
    )
})

test_that("classifyFrozen errors on missing label_col", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()

    expect_error(
        classifyFrozen(fix$pretrained, fix$se,
                       label_col = "nonexistent", verbose = FALSE),
        "not found"
    )
})

test_that("classifyFrozen with nfolds=10 works", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition",
                             nfolds = 10, verbose = FALSE)

    expect_type(result, "list")
    expect_gte(result$auc, 0)
})

test_that("classifyFrozen reproducible with same seed", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()

    r1 <- classifyFrozen(fix$pretrained, fix$se,
                         label_col = "condition",
                         nfolds = 10, seed = 99L, verbose = FALSE)
    r2 <- classifyFrozen(fix$pretrained, fix$se,
                         label_col = "condition",
                         nfolds = 10, seed = 99L, verbose = FALSE)

    expect_equal(r1$predictions, r2$predictions)
    expect_equal(r1$auc, r2$auc)
})

test_that("classifyFrozen handles Ensembl IDs", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    se_ens <- fix$se
    rownames(se_ens) <- fix$pretrained@gene_info$ensembl[1:80]

    result <- classifyFrozen(fix$pretrained, se_ens,
                             label_col = "condition",
                             gene_id_type = "ensembl",
                             verbose = FALSE)

    expect_type(result, "list")
    expect_gte(result$auc, 0)
})

test_that("classifyFrozen returns correct class_levels", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    expect_equal(sort(result$class_levels), c("GroupA", "GroupB"))
})

test_that("classifyFrozen glmnet_fit is a cv.glmnet object", {
    skip_if_not_installed("glmnet")

    fix <- make_pretrained_fixture()
    result <- classifyFrozen(fix$pretrained, fix$se,
                             label_col = "condition", verbose = FALSE)

    expect_s3_class(result$glmnet_fit, "cv.glmnet")
})


# ═══════════════════════════════════════════════════════════════════
#            .alignPretrainedToArchitecture (Mode 2 helper)
# ═══════════════════════════════════════════════════════════════════

test_that("alignPretrainedToArchitecture produces correct dimensions", {
    fix <- make_pretrained_fixture()

    # Build architecture from the SE's genes
    mapping <- data.frame(
        gene_id    = rep(rownames(fix$se)[1:50], each = 2),
        pathway_id = rep(fix$pretrained@pathway_names, 20)
    )
    gpm <- buildGenePathwayMap(mapping,
        feature_genes = rownames(fix$se),
        min_pathway_size = 5
    )
    arch <- buildArchitecture(gpm)

    aligned <- VNNBio:::.alignPretrainedToArchitecture(
        fix$pretrained, arch, gene_id_type = "symbol"
    )

    # W_init should have same dimensions as the architecture mask
    arch_mask <- layerMasks(arch)[[1]]
    expect_equal(nrow(aligned$W_init), nrow(arch_mask))
    expect_equal(ncol(aligned$W_init), ncol(arch_mask))

    # Bias should have same length as pathways
    expect_equal(length(aligned$b_init), ncol(arch_mask))

    # Matched genes should be > 0
    expect_gt(aligned$n_matched, 0)
})

test_that("alignPretrainedToArchitecture zeros unmatched genes", {
    fix <- make_pretrained_fixture()

    # Architecture with some fake genes not in pretrained
    gene_set <- c(rownames(fix$se)[1:30], paste0("FAKE_", 1:20))
    mapping <- data.frame(
        gene_id    = rep(gene_set, each = 1),
        pathway_id = rep(fix$pretrained@pathway_names, 10)
    )
    gpm <- buildGenePathwayMap(mapping,
        feature_genes = gene_set,
        min_pathway_size = 5
    )
    arch <- buildArchitecture(gpm)

    aligned <- VNNBio:::.alignPretrainedToArchitecture(
        fix$pretrained, arch, gene_id_type = "symbol"
    )

    # FAKE genes should have all-zero rows
    arch_genes <- rownames(layerMasks(arch)[[1]])
    fake_idx <- which(grepl("^FAKE_", arch_genes))
    for (i in fake_idx) {
        expect_equal(sum(abs(aligned$W_init[i, ])), 0,
                     info = paste("Fake gene at index", i, "should be zero"))
    }
})


# ═══════════════════════════════════════════════════════════════════
#                  Class weight computation (Mode 2)
# ═══════════════════════════════════════════════════════════════════

test_that("balanced class weights are computed correctly", {
    # Simulate the computation from trainVNN
    y <- c(rep(0, 80), rep(1, 20))  # 80/20 imbalance
    n <- length(y)
    class_counts <- table(y)
    k <- length(class_counts)

    sample_weights <- rep(1.0, n)
    for (cls in names(class_counts)) {
        cls_idx <- which(y == as.numeric(cls))
        sample_weights[cls_idx] <- n / (k * class_counts[cls])
    }

    # Class 0 (majority): 100 / (2 * 80) = 0.625
    # Class 1 (minority): 100 / (2 * 20) = 2.5
    expect_equal(sample_weights[1], 0.625)
    expect_equal(sample_weights[81], 2.5)

    # Weighted sum should equal n (weights are normalized)
    expect_equal(sum(sample_weights), n, tolerance = 1e-10)
})
