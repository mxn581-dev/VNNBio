# =============================================================================
# build_pretrained_weights.R
#
# ONE-TIME script run locally (NOT part of the package).
# Converts ARCHS4-trained numpy weights + R metadata into a single .rds
# bundle that ships in inst/extdata/.
#
# Prereqs:
#   pip install numpy
#   install.packages("reticulate")
#
# Run from the VNNBio_v2 package root:
#   Rscript data-raw/build_pretrained_weights.R
# =============================================================================

library(reticulate)

base_dir  <- "/mnt/Work_SSD_2TB/2026_Harry/archs4/extracted"
pkg_dir   <- "/mnt/Work_SSD_2TB/2026_Harry/VNNBio_v2"

np <- import("numpy")

W_hidden <- np$load(file.path(base_dir, "W_hidden.npy"))
b_hidden <- np$load(file.path(base_dir, "b_hidden.npy"))

gene_info <- readRDS(file.path(base_dir, "frozen_gene_list.rds"))
mask      <- readRDS(file.path(base_dir, "hallmark_mask_aligned.rds"))

# Sanity checks
stopifnot(
    is.matrix(W_hidden) || is.array(W_hidden),
    is.numeric(b_hidden),
    is.data.frame(gene_info),
    all(c("archs4_index", "symbol", "ensembl") %in% colnames(gene_info)),
    nrow(W_hidden) == nrow(gene_info),
    nrow(W_hidden) == nrow(mask),
    ncol(W_hidden) == ncol(mask),
    ncol(W_hidden) == length(b_hidden)
)

# Coerce to standard R types
W_hidden <- as.matrix(W_hidden)
b_hidden <- as.numeric(b_hidden)
mask     <- as.matrix(mask)

# Attach pathway names from mask colnames
pathway_names <- colnames(mask)
if (is.null(pathway_names)) {
    stop("hallmark_mask_aligned.rds must have column names (pathway IDs).")
}

pretrained_bundle <- list(
    W_hidden      = W_hidden,
    b_hidden      = b_hidden,
    mask          = mask,
    gene_info     = gene_info,
    pathway_names = pathway_names,
    metadata      = list(
        weights_id      = "archs4_hallmark_v1",
        source          = "ARCHS4",
        mask_source     = "MSigDB_H (Hallmark)",
        n_genes         = nrow(W_hidden),
        n_pathways      = ncol(W_hidden),
        activation      = "tanh",
        pretrain_epochs = 30L,
        tissue_accuracy = 0.718,
        disease_accuracy = 0.809,
        created         = Sys.time()
    )
)

out_path <- file.path(pkg_dir, "inst", "extdata",
                      "archs4_hallmark_v1.rds")
saveRDS(pretrained_bundle, out_path, compress = "xz")

cat(sprintf(
    "Saved pretrained bundle: %s\n  Genes: %d | Pathways: %d | Size: %.1f MB\n",
    out_path,
    nrow(W_hidden), ncol(W_hidden),
    file.size(out_path) / 1e6
))
