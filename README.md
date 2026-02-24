# VNNBio <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/mxn581-dev/VNNBio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mxn581-dev/VNNBio/actions/workflows/R-CMD-check.yaml)
[![BiocCheck](https://github.com/mxn581-dev/VNNBio/actions/workflows/bioc-check.yaml/badge.svg)](https://github.com/mxn581-dev/VNNBio/actions/workflows/bioc-check.yaml)
[![Codecov](https://codecov.io/gh/mxn581-dev/VNNBio/branch/main/graph/badge.svg)](https://codecov.io/gh/mxn581-dev/VNNBio)
<!-- badges: end -->

**VNNBio** implements Visible Neural Networks (VNNs) for
biologically-constrained interpretable machine learning. Unlike standard
neural networks where hidden layers are opaque, every hidden node in a VNN
corresponds to a known biological concept — a pathway, gene ontology term, or
biological system — making the model's reasoning directly interpretable.

## Why VNNBio?
<img width="676" height="369" alt="vnnbio-removebg-preview" src="https://github.com/user-attachments/assets/cbffd159-3f6c-4a23-8c1f-7fc188920a2c" />

| Feature | VNNBio | DCell | P-NET |
|---|---|---|---|
| **Language** | R + Julia | Python | Python |
| **Bioconductor integration** | ✓ SummarizedExperiment | ✗ | ✗ |
| **Pathway sources** | MSigDB, KEGG, GO, custom | GO only | Reactome |
| **Interpretable output** | Per-pathway importance | Per-GO-term | Per-pathway |

**Key idea:** A sparse mask `M` constrains network connectivity so that gene
*i* can only connect to pathway *j* if gene *i* is annotated to pathway *j*.
After training, pathway importance is read directly from the learned weights —
no post-hoc interpretability methods required.

## Installation

### Prerequisites

VNNBio requires [Julia](https://julialang.org/downloads/) (>= 1.10) for model
training. R-side utilities (mask building, visualization) work without Julia.

```r
# Install from Bioconductor (when available)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VNNBio")

# Or install development version from GitHub
BiocManager::install("mxn581-dev/VNNBio")
```

### Julia backend setup

The first call to `initJuliaBackend()` automatically installs Julia
dependencies (~2-5 minutes on first run).

## Quick Start

```r
library(VNNBio)
library(SummarizedExperiment)

# 1. Load expression data (genes × samples)
se <- your_summarized_experiment

# 2. Build a pathway mask from MSigDB Hallmark gene sets
gpm <- buildMapFromMSigDB(
    category      = "H",
    gene_id_type  = "gene_symbol",
    feature_genes = rownames(se)
)
gpm
#> GenePathwayMap object
#>   Source:     MSigDB_H
#>   Genes:      18,432
#>   Pathways:   50
#>   Nonzeros:   7,321 (0.79% density)

# 3. Define the VNN architecture
arch <- buildArchitecture(gpm, activation = "tanh", n_output = 1L)

# 4. Train
initJuliaBackend()
model <- trainVNN(se, arch, label_col = "condition",
                  task = "classification", epochs = 100L)

# 5. Interpret — which pathways drive the prediction?
importanceScores(model)
#> $`MSigDB_H pathways`
#>     HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION    HALLMARK_ESTROGEN_RESPONSE_EARLY
#>                                          8.432                               6.891 ...

plotImportance(model, top_n = 15)
```

## How It Works
<img width="855" height="52<img width="901" height="657" alt="Screenshot from 2026-02-23 12-28-42" src="https://github.com/user-attachments/assets/0c47f10c-0c55-47f1-a0b9-234dc88711d2" />
6" alt="Screenshot from 2026-02-23 12-23-47" src="https://github.com/user-attachments/assets/8b0a43bb-8d2c-4d20-b45d-fbe70e62413b" />

<img width="900" height="468" alt="Screenshot from 2026-02-23 12-29-30" src="https://github.com/user-attachments/assets/e16f2db3-c9e8-409f-869a-65d41ed52fd9" />
<img width="900" height="468" alt="Screenshot from 2026-02-23 12-29-44" src="https://github.com/user-attachments/assets/d264e4ad-e136-43b0-9474-0d1e732c2108" />


The mask `M` is derived from biological databases (MSigDB, KEGG, GO). Only
connections permitted by the mask can carry signal. After training,
`importanceScores()` reveals which pathways the model relies on.

## Documentation

- **Vignette:** `vignette("VNNBio")` — full walkthrough with real data
- **Reference:** `?trainVNN`, `?buildGenePathwayMap`, `?plotImportance`
- **Architecture guide:** `system.file("ARCHITECTURE.md", package = "VNNBio")`

## Citation

```
@software{VNNBio,
  title  = {VNNBio: Visible Neural Networks for Biologically-Constrained
            Interpretable Machine Learning},
  author = {Minh Hieu Nguyen},
  year   = {2026},
  url    = {https://github.com/mxn581-dev/VNNBio}
}
```

## Related Work

VNNBio builds on the Visible Neural Network paradigm introduced by:

- **DCell** (Ma et al., *Nature Methods* 2018): First VNN mapping the Gene
  Ontology hierarchy to neural network layers.
- **P-NET** (Elmarakeby et al., *Nature* 2021): VNN for predicting prostate
  cancer state from genomic data using Reactome pathways.
- **EMOGI** (Schulte-Sasse et al., *Nature Machine Intelligence* 2021):
  Multi-omics graph integration for cancer gene prediction.

VNNBio is the first R/Bioconductor-native implementation, designed for
integration with standard genomics workflows.

## License

Artistic-2.0
