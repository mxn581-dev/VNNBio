# ================================================================
# VNNBio — Quick Start
#
# Run this script line-by-line in RStudio.
# ================================================================

library(VNNBio)
library(SummarizedExperiment)

# ── Simulated data (swap in your real SE) ──────────────────────
set.seed(42)
n_samples <- 200
n_genes   <- 100

expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
rownames(expr) <- paste0("Gene_", sprintf("%03d", 1:n_genes))
colnames(expr) <- paste0("S", 1:n_samples)

se <- SummarizedExperiment(
    assays  = list(logcounts = expr),
    colData = DataFrame(
        condition = sample(c("Tumor", "Normal"), n_samples, replace = TRUE)
    )
)

# ── Gene-pathway mapping ──────────────────────────────────────
pathway_df <- data.frame(
    gene_id    = c(paste0("Gene_", sprintf("%03d", 1:25)),
                   paste0("Gene_", sprintf("%03d", 15:50)),
                   paste0("Gene_", sprintf("%03d", 40:70)),
                   paste0("Gene_", sprintf("%03d", 60:100))),
    pathway_id = c(rep("Apoptosis", 25), rep("CellCycle", 36),
                   rep("Metabolism", 31), rep("Immune", 41))
)

# ── Build mask (uses accessors, not @ slots) ──────────────────
gpm <- buildGenePathwayMap(pathway_df, feature_genes = rownames(se),
                            min_pathway_size = 5)
gpm  # informative show() output
cat("Density:", round(maskDensity(gpm) * 100, 1), "%\n")

# ── Architecture (convenience constructor) ────────────────────
arch <- buildArchitecture(gpm, activation = "tanh", n_output = 1L)
arch  # shows layer topology

# ── Train (Julia backend) ─────────────────────────────────────
initJuliaBackend()

model <- trainVNN(se, arch,
    label_col     = "condition",
    assay_name    = "logcounts",
    task          = "classification",
    epochs        = 50L,
    learning_rate = 1e-3,
    batch_size    = 32L,
    patience      = 10L      # early stopping
)
model  # informative show() output

# ── Interpret ─────────────────────────────────────────────────
importanceScores(model)          # named list of named vectors
head(importanceTable(model))     # tidy data.frame
plotImportance(model, top_n = 4)
