# =============================================================================
# masks.jl — Sparse mask I/O between R and Julia
# =============================================================================

"""
Global storage for masks and data received from R.
Using a module-level Dict avoids Main namespace pollution.
"""
const MASK_STORE = Dict{String, SparseMatrixCSC{Float32, Int64}}()
const DATA_STORE = Dict{String, Any}()

"""
    receive_sparse_mask(rowval, colptr, nzval, m, n, varname)

Reconstruct a SparseMatrixCSC from the raw CSC components shipped from R's
dgCMatrix. R's `Matrix` package stores:

  - `@i`: 0-based row indices   →  we receive 1-based (converted in R)
  - `@p`: 0-based column pointers →  1-based
  - `@x`: nonzero values

This avoids ever materializing a dense [n_genes × n_pathways] matrix in R
or during transfer — critical when n_genes ~ 20,000 and n_pathways ~ 5,000.
"""
function receive_sparse_mask(
    rowval::Vector{Int64},
    colptr::Vector{Int64},
    nzval::Vector{Float64},
    m::Int64, n::Int64,
    varname::String
)
    mask = SparseMatrixCSC{Float32, Int64}(
        m, n,
        colptr,
        rowval,
        Float32.(nzval)
    )
    MASK_STORE[varname] = mask
    @info "Received mask '$varname': $(m)×$(n), nnz=$(nnz(mask)), " *
          "density=$(round(nnz(mask) / (m * n) * 100; digits=2))%"
    return nothing
end


"""
    receive_data(X, y, n_output, activation)

Receive the expression matrix and labels from R.

- `X`: Matrix of shape [n_samples × n_genes] (already transposed in R).
- `y`: Vector of labels.
- `n_output`: Number of output neurons.
- `activation`: Activation function name ("tanh", "relu", "sigmoid").
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
