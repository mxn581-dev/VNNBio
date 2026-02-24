# VNNBio 0.99.0

## New Features

* Initial Bioconductor submission.
* `GenePathwayMap` S4 class for encoding gene-pathway associations as sparse
  adjacency matrices.
* `VNNArchitecture` S4 class for defining layer-wise network topology with
  dimensional chain validation.
* `VNNModel` S4 class for storing trained models with provenance tracking.
* `buildGenePathwayMap()` constructs sparse masks from gene-pathway
  data.frames with automatic ID deduplication, pathway size filtering, and
  feature alignment to SummarizedExperiment objects.
* `buildMapFromMSigDB()` fetches gene sets directly from MSigDB via the
  msigdbr package.
* `harmonizeGeneIds()` translates between Ensembl, Entrez, and Symbol ID
  types via org.Hs.eg.db.
* `trainVNN()` dispatches model training to a Julia/Lux.jl backend via
  JuliaConnectoR with biology-aware L1 regularization.
* `predict()` method for `VNNModel` objects.
* `plotImportance()`, `plotWeightHeatmap()`, and `plotArchitecture()` for
  publication-ready interpretability visualization.
* Full S4 accessor functions for all classes.
* `buildArchitecture()` convenience method for single-layer architectures.
* `importanceTable()` returns pathway scores as a tidy data.frame.
