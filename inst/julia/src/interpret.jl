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
