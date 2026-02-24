# =============================================================================
# model.jl — Custom Lux layers for biologically-constrained connectivity
# =============================================================================

"""
    MaskedDense(mask::AbstractMatrix; activation=tanh)

A Lux-compatible dense layer where the weight matrix is elementwise-multiplied
by a fixed binary mask at every forward pass. This enforces the biological
constraint: gene i can only influence pathway j if mask[i, j] == 1.

## Why Hadamard masking (W .* M) instead of sparse storage?

1. **Gradient flow**: Autodiff (Zygote) works seamlessly on dense W .* M.
   Directly parameterizing only the nonzero positions would require custom
   adjoints and complicates optimizer state management.

2. **GPU compatibility**: Dense GEMM with a mask is far better optimized on
   GPUs than truly sparse matmuls for the density levels typical in biology
   (5–30% nonzero).

3. **Interpretability**: After training, we can inspect the full W .* M to
   see both active and suppressed connections.
"""
struct MaskedDense{M <: AbstractMatrix, F} <: Lux.AbstractLuxLayer
    mask::M
    in_dims::Int
    out_dims::Int
    activation::F
end

function MaskedDense(
    mask::AbstractMatrix;
    activation::Function = tanh
)
    in_dims, out_dims = size(mask)
    return MaskedDense(
        Float32.(Matrix(mask)),  # densify for GPU-friendly ops
        in_dims,
        out_dims,
        activation
    )
end

"""
    Lux.initialparameters(rng, layer::MaskedDense)

Xavier-uniform initialization, pre-masked so that zeroed connections start
at zero and the optimizer never wastes capacity on them.
"""
function Lux.initialparameters(rng::AbstractRNG, layer::MaskedDense)
    # Xavier uniform scaled by mask density for stable initialization
    scale = sqrt(2.0f0 / (layer.in_dims + layer.out_dims))
    W = randn(rng, Float32, layer.in_dims, layer.out_dims) .* scale
    W = W .* layer.mask  # zero out forbidden connections from the start
    b = zeros(Float32, layer.out_dims)
    return (weight = W, bias = b)
end

function Lux.initialstates(::AbstractRNG, layer::MaskedDense)
    return (mask = layer.mask,)
end

function Lux.parameterlength(layer::MaskedDense)
    return layer.in_dims * layer.out_dims + layer.out_dims
end

"""
    (layer::MaskedDense)(x, ps, st)

Forward pass: z = activation.( (W .* mask)' * x .+ b )

Input `x` has shape [in_dims, batch_size].
Output has shape [out_dims, batch_size].
"""
function (layer::MaskedDense)(x::AbstractMatrix, ps, st)
    W_masked = ps.weight .* st.mask   # enforce sparsity every forward pass
    z = layer.activation.(W_masked' * x .+ ps.bias)
    return z, st
end


# =============================================================================
# Model builder: stack MaskedDense layers + final Dense output
# =============================================================================

"""
    build_vnn_model(mask_names; n_output=1, activation_fn=tanh)

Construct a Lux Chain from an ordered list of mask names (keys into
MASK_STORE). The final layer is a standard Dense layer to the output.

Returns a `Lux.Chain` model.
"""
function build_vnn_model(
    mask_names::Vector{String};
    n_output::Int = 1,
    activation_fn::Function = tanh
)
    layers = Lux.AbstractLuxLayer[]

    for (i, name) in enumerate(mask_names)
        mask = MASK_STORE[name]
        push!(layers, MaskedDense(mask; activation = activation_fn))
    end

    # Final dense layer from last hidden size to output
    last_mask = MASK_STORE[mask_names[end]]
    last_hidden = size(last_mask, 2)

    if n_output == 1
        # Binary classification or scalar regression
        push!(layers, Lux.Dense(last_hidden, 1, sigmoid))
    else
        push!(layers, Lux.Dense(last_hidden, n_output))
    end

    return Lux.Chain(layers...)
end


"""
    get_activation_fn(name::String) → Function

Maps activation function name strings from R to Julia functions.
"""
function get_activation_fn(name::String)
    activations = Dict(
        "tanh"    => tanh,
        "relu"    => relu,
        "sigmoid" => sigmoid,
        "gelu"    => gelu,
        "swish"   => swish,
    )
    return get(activations, lowercase(name)) do
        @warn "Unknown activation '$name', defaulting to tanh"
        tanh
    end
end
