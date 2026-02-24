# =============================================================================
# visualization.R -- Publication-ready plots for VNN interpretability
# =============================================================================

# =============================================================================
# Bar plot of top-N pathway importance per layer
# =============================================================================

#' Plot Pathway Importance Scores
#'
#' Horizontal bar chart of the top-N most important pathways (or nodes) at
#' each biological layer of the VNN.
#'
#' @param model A trained \code{VNNModel} object.
#' @param top_n Integer. Number of top pathways to show per layer. Default 20.
#' @param normalize Logical. Scale scores to the range 0 to 1 within each
#'   layer? Default TRUE.
#'
#' @return A \code{ggplot} object.
#' @export
plotImportance <- function(model, top_n = 20L, normalize = TRUE) {

    scores <- importanceScores(model)
    if (length(scores) == 0) {
        stop("No importance scores found. Was the model trained?")
    }

    plot_df <- do.call(rbind, lapply(names(scores), function(layer) {
        s <- scores[[layer]]
        if (normalize && max(s) > 0) s <- s / max(s)
        top_idx <- seq_len(min(top_n, length(s)))
        data.frame(
            layer   = layer,
            pathway = names(s)[top_idx],
            score   = unname(s[top_idx]),
            stringsAsFactors = FALSE
        )
    }))

    plot_df$pathway <- factor(plot_df$pathway,
                              levels = rev(unique(plot_df$pathway)))

    ggplot(plot_df, aes(x = .data$pathway, y = .data$score,
                        fill = .data$score)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7",
                             high = "#b2182b", midpoint = 0.5) +
        facet_wrap(~ layer, scales = "free_y") +
        labs(
            title = "VNN Pathway Importance Scores",
            x = NULL,
            y = "Importance (normalized sum |W * mask|)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            strip.text = element_text(face = "bold"),
            axis.text.y = element_text(size = 8)
        )
}


# =============================================================================
# Heatmap of the masked weight matrix
# =============================================================================

#' Plot Masked Weight Heatmap
#'
#' Visualizes the actual trained weight values at a specific layer, masked by
#' the biological connectivity. Uses ComplexHeatmap for publication quality.
#'
#' @param model A trained \code{VNNModel} object.
#' @param layer Integer. Which layer's weights to display (1-indexed).
#' @param top_n_rows Integer. Show only top N input features by total
#'   outgoing weight magnitude. Default 50.
#' @param top_n_cols Integer. Show only top N output nodes by incoming weight
#'   magnitude. Default 30.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object (draws automatically).
#' @export
plotWeightHeatmap <- function(model, layer = 1L,
                              top_n_rows = 50L, top_n_cols = 30L) {

    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("Install ComplexHeatmap:  BiocManager::install('ComplexHeatmap')")
    }
    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Install circlize:  install.packages('circlize')")
    }

    ## Pull weight matrix from Julia using model_ref for correct routing
    model_ref <- model@julia_model_ref
    W <- as.matrix(
        JuliaConnectoR::juliaCall("_vnn_get_masked_weights",
                                  as.integer(layer), model_ref)
    )

    ## Subset to most important rows/columns
    row_imp <- rowSums(abs(W))
    col_imp <- colSums(abs(W))

    keep_rows <- order(row_imp, decreasing = TRUE)[
        seq_len(min(top_n_rows, nrow(W)))]
    keep_cols <- order(col_imp, decreasing = TRUE)[
        seq_len(min(top_n_cols, ncol(W)))]

    W_sub <- W[keep_rows, keep_cols]

    ## Apply mask dimnames if available
    arch <- architecture(model)
    masks <- layerMasks(arch)
    if (layer <= length(masks)) {
        mask <- masks[[layer]]
        if (!is.null(rownames(mask)))
            rownames(W_sub) <- rownames(mask)[keep_rows]
        if (!is.null(colnames(mask)))
            colnames(W_sub) <- colnames(mask)[keep_cols]
    }

    col_fun <- circlize::colorRamp2(
        c(-max(abs(W_sub)), 0, max(abs(W_sub))),
        c("#2166ac", "#f7f7f7", "#b2182b")
    )

    ComplexHeatmap::Heatmap(
        W_sub,
        name = "Weight",
        col = col_fun,
        column_title = paste0("Layer ", layer, " - Masked Weights"),
        row_title = "Input Features",
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        row_names_gp = grid::gpar(fontsize = 7),
        column_names_gp = grid::gpar(fontsize = 8),
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2"
    )
}


# =============================================================================
# Network topology visualization with igraph
# =============================================================================

#' Plot VNN Architecture as a Layered Graph
#'
#' Renders the VNN topology as a bipartite/multipartite graph where node
#' size represents importance and edge width represents |weight|.
#'
#' @param model A trained \code{VNNModel} object.
#' @param top_n_per_layer Integer. Max nodes shown per layer. Default 15.
#' @param edge_threshold Numeric. Only show edges with |weight| above this.
#'   Default 0.01.
#'
#' @return An \code{igraph} plot (drawn as side effect).
#' @export
plotArchitecture <- function(model, top_n_per_layer = 15L,
                             edge_threshold = 0.01) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Install igraph:  install.packages('igraph')")
    }
    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("Install scales:  install.packages('scales')")
    }

    arch <- architecture(model)
    scores <- importanceScores(model)
    masks <- layerMasks(arch)
    lnames <- layerNames(arch)

    nodes <- data.frame(name = character(), layer = integer(),
                        importance = numeric(), stringsAsFactors = FALSE)
    edges <- data.frame(from = character(), to = character(),
                        weight = numeric(), stringsAsFactors = FALSE)

    for (i in seq_along(masks)) {
        mask <- masks[[i]]
        layer_label <- if (i <= length(lnames)) lnames[i] else paste0("layer_", i)

        ## Determine top output nodes for this layer
        if (layer_label %in% names(scores)) {
            imp <- scores[[layer_label]]
            top_idx <- seq_len(min(top_n_per_layer, length(imp)))
            top_out <- names(imp)[top_idx]
            out_imp <- imp[top_out]
        } else {
            top_out <- colnames(mask)[seq_len(min(top_n_per_layer, ncol(mask)))]
            out_imp <- setNames(rep(1, length(top_out)), top_out)
        }

        ## Tag output node names with layer to avoid name collisions
        out_tagged <- paste0(layer_label, "::", top_out)

        ## Determine input nodes connected to these outputs
        ## For first layer: inputs are genes; for later layers: previous outputs
        if (i == 1) {
            in_prefix <- "input"
        } else {
            prev_label <- if (i - 1 <= length(lnames)) lnames[i - 1] else paste0("layer_", i - 1)
            in_prefix <- prev_label
        }

        ## Build edges from mask structure
        ## mask is [n_input x n_output]: mask[gene, pathway] = 1
        top_out_idx <- match(top_out, colnames(mask))

        for (ci in seq_along(top_out_idx)) {
            col_j <- top_out_idx[ci]
            if (is.na(col_j)) next

            ## Find nonzero input rows for this output node
            col_start <- mask@p[col_j] + 1L
            col_end <- mask@p[col_j + 1L]
            if (col_end < col_start) next

            row_indices <- mask@i[col_start:col_end] + 1L  # 0-based to 1-based
            input_names <- rownames(mask)[row_indices]

            ## Subsample inputs if too many
            if (length(input_names) > top_n_per_layer) {
                input_names <- sample(input_names, top_n_per_layer)
            }

            in_tagged <- paste0(in_prefix, "::", input_names)

            ## Add input nodes (if first layer)
            if (i == 1) {
                new_inputs <- setdiff(in_tagged, nodes$name)
                if (length(new_inputs) > 0) {
                    nodes <- rbind(nodes, data.frame(
                        name = new_inputs, layer = 0L,
                        importance = 0.5,
                        stringsAsFactors = FALSE))
                }
            }

            ## Add edges
            new_edges <- data.frame(
                from   = in_tagged,
                to     = rep(out_tagged[ci], length(in_tagged)),
                weight = 1.0,
                stringsAsFactors = FALSE
            )
            edges <- rbind(edges, new_edges)
        }

        ## Add output nodes
        nodes <- rbind(nodes, data.frame(
            name       = out_tagged,
            layer      = i,
            importance = unname(out_imp),
            stringsAsFactors = FALSE
        ))
    }

    ## Add output layer node
    nodes <- rbind(nodes, data.frame(
        name = "Output", layer = length(masks) + 1L,
        importance = 1.0, stringsAsFactors = FALSE))

    ## Connect last hidden layer to output
    last_outputs <- nodes$name[nodes$layer == length(masks)]
    edges <- rbind(edges, data.frame(
        from = last_outputs, to = rep("Output", length(last_outputs)),
        weight = 1.0, stringsAsFactors = FALSE))

    ## Remove duplicate nodes
    nodes <- nodes[!duplicated(nodes$name), ]

    ## Filter edges by threshold
    edges <- edges[edges$weight >= edge_threshold, ]

    ## Build graph
    g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

    ## Style
    n_layers <- max(nodes$layer) + 1
    layer_colors <- grDevices::hcl.colors(n_layers, "Set2")
    igraph::V(g)$color <- layer_colors[nodes$layer + 1]
    igraph::V(g)$size <- scales::rescale(nodes$importance, to = c(3, 12))

    ## Strip prefix for labels (show only the gene/pathway name)
    clean_labels <- sub("^.*::", "", igraph::V(g)$name)
    igraph::V(g)$label <- clean_labels
    igraph::V(g)$label.cex <- 0.5

    ## Layered layout
    layout <- igraph::layout_with_sugiyama(g)$layout

    plot(g, layout = layout, edge.arrow.size = 0.2,
         edge.color = grDevices::adjustcolor("gray40", alpha.f = 0.5),
         main = "VNN Architecture")

    invisible(g)
}
