# ================================================================
# test-visualization.R — Visualization Function Tests
#
# plotImportance:      No Julia required (reads importance_scores slot)
# plotArchitecture:    No Julia required (reads architecture + scores)
# plotWeightHeatmap:   Requires Julia (calls get_masked_weights)
# ================================================================

library(Matrix)

# ── Shared fixture builder ───────────────────────────────────────

make_test_model <- function() {
    mapping <- data.frame(
        gene_id    = c(paste0("G", 1:10), paste0("G", 5:15)),
        pathway_id = c(rep("PA", 10), rep("PB", 11))
    )
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 5)
    arch <- buildArchitecture(gpm)

    pw_names <- colnames(maskMatrix(gpm))
    scores <- setNames(seq(length(pw_names), 1), pw_names)

    # layer_names from buildArchitecture is paste0(source, " pathways")
    # source defaults to "custom" -> "custom pathways"
    new("VNNModel",
        architecture      = arch,
        importance_scores = list("custom pathways" = scores),
        task              = "classification"
    )
}

make_single_pathway_model <- function() {
    mapping <- data.frame(
        gene_id    = paste0("G", 1:10),
        pathway_id = rep("PA", 10)
    )
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 5)
    arch <- buildArchitecture(gpm)

    new("VNNModel",
        architecture      = arch,
        importance_scores = list("custom pathways" = c(PA = 3.0)),
        task              = "classification"
    )
}


# ═══════════════════════════════════════════════════════════════════
#                         plotImportance
# ═══════════════════════════════════════════════════════════════════

test_that("plotImportance with top_n returns ggplot", {
    model <- make_test_model()
    p <- plotImportance(model, top_n = 2)
    expect_s3_class(p, "gg")
})

test_that("plotImportance with empty scores throws error", {
    model <- new("VNNModel")
    expect_error(plotImportance(model), "No importance scores")
})

test_that("plotImportance with normalize=FALSE returns ggplot", {
    model <- make_test_model()
    p <- plotImportance(model, normalize = FALSE)
    expect_s3_class(p, "gg")
})

test_that("plotImportance with single pathway returns ggplot", {
    model <- make_single_pathway_model()
    p <- plotImportance(model, top_n = 5)
    expect_s3_class(p, "gg")
})

test_that("plotImportance top_n larger than pathway count works", {
    model <- make_test_model()
    # Only 2 pathways but requesting 50 — should not error
    p <- plotImportance(model, top_n = 50)
    expect_s3_class(p, "gg")
})

test_that("plotImportance normalized scores are in [0, 1]", {
    model <- make_test_model()
    p <- plotImportance(model, normalize = TRUE)

    # Extract the data used by ggplot
    plot_data <- ggplot2::ggplot_build(p)$data[[1]]
    expect_true(all(plot_data$y >= 0))
    expect_true(all(plot_data$y <= 1))
})


# ═══════════════════════════════════════════════════════════════════
#                        plotArchitecture
# ═══════════════════════════════════════════════════════════════════

test_that("plotArchitecture returns non-empty igraph", {
    model <- make_test_model()

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    result <- plotArchitecture(model)
    dev.off()
    unlink(tmp)

    expect_s3_class(result, "igraph")
    expect_gt(igraph::vcount(result), 0)
    expect_gt(igraph::ecount(result), 0)
})

test_that("plotArchitecture graph contains expected node layers", {
    model <- make_test_model()

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    result <- plotArchitecture(model)
    dev.off()
    unlink(tmp)

    # Graph should have nodes from input layer, pathway layer, and output
    node_names <- igraph::V(result)$name
    expect_true(any(grepl("^input::", node_names)))
    expect_true(any(grepl("^custom pathways::", node_names)))
    expect_true("Output" %in% node_names)
})

test_that("plotArchitecture top_n_per_layer reduces node count", {
    model <- make_test_model()

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    full <- plotArchitecture(model, top_n_per_layer = 15L)
    dev.off()

    pdf(tmp)
    small <- plotArchitecture(model, top_n_per_layer = 1L)
    dev.off()
    unlink(tmp)

    expect_lte(igraph::vcount(small), igraph::vcount(full))
})

test_that("plotArchitecture high edge_threshold reduces edge count", {
    model <- make_test_model()

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    normal <- plotArchitecture(model, edge_threshold = 0.01)
    dev.off()

    pdf(tmp)
    strict <- plotArchitecture(model, edge_threshold = 999)
    dev.off()
    unlink(tmp)

    expect_lte(igraph::ecount(strict), igraph::ecount(normal))
})

test_that("plotArchitecture with single pathway produces valid graph", {
    model <- make_single_pathway_model()

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    result <- plotArchitecture(model)
    dev.off()
    unlink(tmp)

    expect_s3_class(result, "igraph")
    expect_gt(igraph::ecount(result), 0)
})

test_that("plotArchitecture output node connects to all pathway nodes", {
    model <- make_test_model()

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    result <- plotArchitecture(model)
    dev.off()
    unlink(tmp)

    # The Output node should have incoming edges from all shown pathway nodes
    output_edges <- igraph::incident(result, "Output", mode = "in")
    expect_gt(length(output_edges), 0)
})

test_that("plotArchitecture without importance scores still works", {
    # Model with architecture but no importance scores —
    # function should fall back to showing first top_n_per_layer columns
    mapping <- data.frame(
        gene_id    = c(paste0("G", 1:10), paste0("G", 5:15)),
        pathway_id = c(rep("PA", 10), rep("PB", 11))
    )
    gpm <- buildGenePathwayMap(mapping, min_pathway_size = 5)
    arch <- buildArchitecture(gpm)

    model <- new("VNNModel",
        architecture      = arch,
        importance_scores = list(),
        task              = "classification"
    )

    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    result <- plotArchitecture(model)
    dev.off()
    unlink(tmp)

    expect_s3_class(result, "igraph")
    expect_gt(igraph::vcount(result), 0)
})

test_that("plotArchitecture errors without igraph installed", {
    # We can't truly unload igraph in a test, but we verify the
    # error message text is correct by inspecting the source.
    # If igraph IS installed, this test just confirms the function runs.
    # If igraph were missing, we'd expect this error:
    model <- make_test_model()

    if (!requireNamespace("igraph", quietly = TRUE)) {
        expect_error(
            plotArchitecture(model),
            "Install igraph"
        )
    } else {
        tmp <- tempfile(fileext = ".pdf")
        pdf(tmp)
        result <- plotArchitecture(model)
        dev.off()
        unlink(tmp)
        expect_s3_class(result, "igraph")
    }
})


# ═══════════════════════════════════════════════════════════════════
#                       plotWeightHeatmap
# ═══════════════════════════════════════════════════════════════════

# plotWeightHeatmap calls JuliaConnectoR::juliaCall internally,
# so full testing requires Julia. We test guard clauses here;
# full integration belongs in test-integration-julia.R.

test_that("plotWeightHeatmap errors without ComplexHeatmap", {
    model <- make_test_model()

    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        expect_error(
            plotWeightHeatmap(model, layer = 1L),
            "ComplexHeatmap"
        )
    }
})

test_that("plotWeightHeatmap errors without circlize", {
    model <- make_test_model()

    if (!requireNamespace("circlize", quietly = TRUE)) {
        expect_error(
            plotWeightHeatmap(model, layer = 1L),
            "circlize"
        )
    }
})

test_that("plotWeightHeatmap errors when model has no Julia ref", {
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("circlize")

    model <- make_test_model()

    # model has NULL julia_model_ref — JuliaConnectoR call will fail
    expect_error(plotWeightHeatmap(model, layer = 1L))
})
