# =============================================================================
# masks.jl — Sparse mask I/O between R and Julia
#
# PATCHED: Added receive_init_weights() for Mode 2 pretrained initialization.
# Replace the existing masks.jl with this file.
# =============================================================================

"""
Global storage for masks and data received from R.
"""
const MASK_STORE = Dict{String, SparseMatrixCSC{Float32, Int64}}()
const DATA_STORE = Dict{String, Any}()
const INIT_STORE = Dict{String, Any}()  # NEW: pretrained weight storage

"""
    receive_sparse_mask(rowval, colptr, nzval, m, n, varname)

Reconstruct a SparseMatrixCSC from the raw CSC components shipped from R's
dgCMatrix.
"""
function receive_sparse_mask(
    rowval::Vector{Int64},
    colptr::Vector{Int64},
    nzval::Vector{Float64},
    m::Int64, n::Int64,
    varname::String
)
    mask = SparseMatrixCSC{Float32, Int64}(
        m, n, colptr, rowval, Float32.(nzval)
    )
    MASK_STORE[varname] = mask
    @info "Received mask '$varname': $(m)×$(n), nnz=$(nnz(mask)), " *
          "density=$(round(nnz(mask) / (m * n) * 100; digits=2))%"
    return nothing
end


"""
    receive_data(X, y, n_output, activation)

Receive the expression matrix and labels from R.
"""
function receive_data(
    X::Matrix{Float64},
    y::Vector{Float64},
    n_output::Int64,
    activation::String
)
    DATA_STORE["X"] = Float32.(X)
    DATA_STORE["y"] = Float32.(y)
    DATA_STORE["n_output"] = n_output
    DATA_STORE["activation"] = activation
    @info "Received data: $(size(X, 1)) samples × $(size(X, 2)) features, " *
          "$(n_output) output(s), activation=$(activation)"
    return nothing
end


"""
    receive_init_weights(W_init, b_init)

Receive pretrained weights from R for Mode 2 (fine-tuned) initialization.
These replace the random Xavier init for the first hidden layer.

- `W_init`: Matrix [n_genes × n_pathways] of pretrained weights.
- `b_init`: Vector [n_pathways] of pretrained biases.
"""
function receive_init_weights(
    W_init::Matrix{Float64},
    b_init::Vector{Float64}
)
    INIT_STORE["W_init"] = Float32.(W_init)
    INIT_STORE["b_init"] = Float32.(b_init)
    @info "Received pretrained init weights: $(size(W_init, 1))×$(size(W_init, 2))"
    return nothing
end
