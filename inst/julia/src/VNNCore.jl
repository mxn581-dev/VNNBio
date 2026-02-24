"""
    VNNCore

Julia backend for the VNNBio R package. Implements biologically-constrained
Visible Neural Networks using Lux.jl for the model definition and Optimisers.jl
for training. All masks are received from R as sparse CSC components and
assembled into native Julia SparseMatrixCSC objects.

Architecture:

    Input (genes) ──[mask₁]──▶ Hidden Layer 1 (pathways)
                               ──[mask₂]──▶ Hidden Layer 2 (systems)
                                             ──[dense]──▶ Output

Each masked layer performs:  z = activation.( (W .* mask) * x .+ b )
where .* is Hadamard (elementwise) product enforcing sparsity.
"""
module VNNCore

using SparseArrays
using Statistics
using Random
using Lux
using Optimisers
using Zygote

export receive_sparse_mask, receive_data, train_vnn, get_masked_weights, predict_vnn

include("masks.jl")
include("model.jl")
include("training.jl")
include("interpret.jl")

end # module VNNCore
