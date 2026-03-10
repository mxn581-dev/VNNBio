# =============================================================================
# interpret.jl — Weight extraction and prediction with model_ref routing
# =============================================================================

"""
    get_masked_weights(layer_index, model_ref)

Extract W .* M for the specified layer from a specific trained model.
If model_ref is empty or missing, falls back to the last trained model.
"""
function get_masked_weights(layer_index::Int64, model_ref::String="")
    if isempty(MODEL_STORE)
        error("No trained model found. Call train_vnn first.")
    end

    model_key = if isempty(model_ref)
        last(collect(keys(MODEL_STORE)))
    else
        if !haskey(MODEL_STORE, model_ref)
            error("Model '$model_ref' not found in MODEL_STORE. " *
                  "Available: $(collect(keys(MODEL_STORE)))")
        end
        model_ref
    end

    stored = MODEL_STORE[model_key]
    ps = stored["ps"]
    st = stored["st"]

    layer_keys = collect(keys(ps))
    if layer_index > length(layer_keys)
        error("Layer index $layer_index exceeds number of layers ($(length(layer_keys)))")
    end

    lk = layer_keys[layer_index]
    W = ps[lk].weight

    if hasproperty(st[lk], :mask)
        return Matrix(W .* st[lk].mask)
    else
        return Matrix(W)
    end
end


"""
    extract_vnn_params(model_ref)

Extract all model parameters (weights, biases, masks) for R-side computation.
Returns a Dict with:
  - "n_layers": total number of layers
  - "layer_<i>_weight": weight matrix for layer i
  - "layer_<i>_bias": bias vector for layer i
  - "layer_<i>_has_mask": whether layer i has a mask (MaskedDense vs Dense)
  - "layer_<i>_mask": mask matrix (only if has_mask is true)

This enables pure-R forward passes for Shapley attribution without
repeated Julia calls per coalition evaluation.
"""
function extract_vnn_params(model_ref::String="")
    if isempty(MODEL_STORE)
        error("No trained model found. Call train_vnn first.")
    end

    model_key = if isempty(model_ref)
        last(collect(keys(MODEL_STORE)))
    else
        if !haskey(MODEL_STORE, model_ref)
            error("Model '$model_ref' not found in MODEL_STORE. " *
                  "Available: $(collect(keys(MODEL_STORE)))")
        end
        model_ref
    end

    stored = MODEL_STORE[model_key]
    ps = stored["ps"]
    st = stored["st"]

    layer_keys = collect(keys(ps))
    n_layers = length(layer_keys)

    result = Dict{String, Any}()
    result["n_layers"] = n_layers

    for (i, lk) in enumerate(layer_keys)
        result["layer_$(i)_weight"] = Matrix(Float64.(ps[lk].weight))
        result["layer_$(i)_bias"]   = Vector(Float64.(ps[lk].bias))

        has_mask = hasproperty(st[lk], :mask)
        result["layer_$(i)_has_mask"] = has_mask

        if has_mask
            result["layer_$(i)_mask"] = Matrix(Float64.(st[lk].mask))
        end
    end

    @info "Extracted parameters for $n_layers layers from model '$model_key'"
    return result
end


"""
    predict_vnn(X_new, model_ref)

Generate predictions from a specific trained model.
"""
function predict_vnn(X_new::Matrix{Float64}, model_ref::String="")
    if isempty(MODEL_STORE)
        error("No trained model found.")
    end

    model_key = if isempty(model_ref)
        last(collect(keys(MODEL_STORE)))
    else
        if !haskey(MODEL_STORE, model_ref)
            error("Model '$model_ref' not found in MODEL_STORE.")
        end
        model_ref
    end

    stored = MODEL_STORE[model_key]
    y_hat, _ = stored["model"](Float32.(X_new'), stored["ps"],
                                Lux.testmode(stored["st"]))
    return Matrix(y_hat)
end
