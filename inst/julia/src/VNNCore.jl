"""
    VNNCore

Julia backend for the VNNBio R package. Implements biologically-constrained
Visible Neural Networks using Lux.jl for the model definition and Optimisers.jl
for training.

PATCHED: Added receive_init_weights export for Mode 2 pretrained init.
"""
module VNNCore

using SparseArrays
using Statistics
using Random
using Lux
using Optimisers
using Zygote

export receive_sparse_mask, receive_data, receive_init_weights,
       train_vnn, get_masked_weights, predict_vnn, extract_vnn_params

include("masks.jl")
include("model.jl")
include("training.jl")
include("interpret.jl")

end # module VNNCore
