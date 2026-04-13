# =============================================================================
# training.jl — Training loop with biology-aware regularization + early stopping
#
# PATCHED vs original:
#   1. train_vnn() accepts sample_weights vector for class-weighted loss
#   2. train_vnn() accepts use_pretrained flag to init from INIT_STORE
#   3. Backward-compatible: old 10-arg call signature still works
#
# Replace the existing training.jl with this file.
# =============================================================================

const MODEL_STORE = Dict{String, Any}()

# Backward-compatible: old call with 10 args (no sample_weights, no pretrained)
function train_vnn(mask_name::String, args...)
    return train_vnn([mask_name], args...)
end

# Old 10-arg signature: inject default sample_weights and use_pretrained=false
function train_vnn(
    mask_names::Vector{String},
    epochs::Int64,
    lr::Float64,
    batch_size::Int64,
    val_frac::Float64,
    l1_lambda::Float64,
    patience::Int64,
    min_delta::Float64,
    seed::Int64,
    verbose::Bool
)
    n_samples = size(DATA_STORE["X"], 1)
    sample_weights = ones(Float64, n_samples)
    return train_vnn(
        mask_names, epochs, lr, batch_size, val_frac, l1_lambda,
        patience, min_delta, seed, verbose, sample_weights, false
    )
end

# Full 12-arg signature with class weighting + pretrained init
function train_vnn(
    mask_names::Vector{String},
    epochs::Int64,
    lr::Float64,
    batch_size::Int64,
    val_frac::Float64,
    l1_lambda::Float64,
    patience::Int64,
    min_delta::Float64,
    seed::Int64,
    verbose::Bool,
    sample_weights::Vector{Float64},
    use_pretrained::Bool
)
    rng = Random.MersenneTwister(seed)

    X_all = DATA_STORE["X"]
    y_all = DATA_STORE["y"]
    n_output = DATA_STORE["n_output"]
    act_name = DATA_STORE["activation"]
    sw_all = Float32.(sample_weights)

    n_samples = size(X_all, 1)

    # Train/val split
    n_val = max(1, round(Int, n_samples * val_frac))
    perm = randperm(rng, n_samples)
    val_idx = perm[1:n_val]
    train_idx = perm[n_val+1:end]

    X_train = X_all[train_idx, :]'
    y_train = reshape(y_all[train_idx], 1, :)
    sw_train = sw_all[train_idx]
    X_val   = X_all[val_idx, :]'
    y_val   = reshape(y_all[val_idx], 1, :)
    sw_val  = sw_all[val_idx]
    n_train = length(train_idx)

    # Build model
    act_fn = get_activation_fn(act_name)
    model = build_vnn_model(mask_names; n_output = n_output, activation_fn = act_fn)
    ps, st = Lux.setup(rng, model)

    # ── Pretrained initialization (Mode 2) ──
    if use_pretrained && haskey(INIT_STORE, "W_init")
        W_init = INIT_STORE["W_init"]
        b_init = INIT_STORE["b_init"]

        layer_keys = collect(keys(ps))
        first_key = layer_keys[1]
        orig_W = ps[first_key].weight
        orig_b = ps[first_key].bias

        # Dimension check: W_init must match first layer weight shape
        if size(W_init) == size(orig_W)
            # Apply mask to pretrained weights (zero out forbidden connections)
            if hasproperty(st[first_key], :mask)
                W_init_masked = W_init .* st[first_key].mask
            else
                W_init_masked = W_init
            end

            # Replace the first layer's parameters
            new_first = (weight = W_init_masked, bias = Float32.(b_init))
            ps = merge(ps, NamedTuple{(first_key,)}((new_first,)))

            if verbose
                @info "Initialized first layer from pretrained weights " *
                      "($(size(W_init, 1))×$(size(W_init, 2)))"
            end
        else
            @warn "Pretrained weight dimensions $(size(W_init)) don't match " *
                  "layer dimensions $(size(orig_W)). Using random init."
        end
    end

    # Optimizer
    opt = Optimisers.Adam(Float32(lr))
    opt_state = Optimisers.setup(opt, ps)

    # Config
    l1_weight = Float32(l1_lambda)
    is_cls = (n_output == 1)
    use_early_stop = patience > 0

    train_losses = Float32[]
    val_losses   = Float32[]

    # Early stopping state
    best_val_loss = Inf32
    best_ps = ps
    wait = 0

    for epoch in 1:epochs
        shuf = randperm(rng, n_train)
        epoch_loss = 0.0f0
        n_batches = 0

        for bs in 1:batch_size:n_train
            be = min(bs + batch_size - 1, n_train)
            idx = shuf[bs:be]
            xb = X_train[:, idx]
            yb = y_train[:, idx]
            wb = sw_train[idx]  # per-sample weights for this batch

            result = Zygote.withgradient(ps) do p
                y_hat, _ = model(xb, p, st)

                if is_cls
                    eps = 1.0f-7
                    yc = clamp.(y_hat, eps, 1.0f0 - eps)
                    # Weighted binary cross-entropy
                    per_sample_loss = -(yb .* log.(yc) .+ (1.0f0 .- yb) .* log.(1.0f0 .- yc))
                    # Apply sample weights and average
                    data_loss = sum(per_sample_loss .* reshape(wb, 1, :)) / sum(wb)
                else
                    per_sample_loss = (y_hat .- yb) .^ 2
                    data_loss = sum(per_sample_loss .* reshape(wb, 1, :)) / sum(wb)
                end

                reg = 0.0f0
                for lp in values(p)
                    if lp isa NamedTuple && haskey(lp, :weight)
                        reg = reg + sum(abs.(lp.weight))
                    end
                end

                data_loss + l1_weight * reg
            end

            loss_val = result.val
            grad_ps = result.grad[1]

            opt_state, ps = Optimisers.update(opt_state, ps, grad_ps)
            epoch_loss += loss_val
            n_batches += 1
        end

        push!(train_losses, epoch_loss / n_batches)

        # Validation loss (also weighted for consistency)
        y_hat_val, _ = model(X_val, ps, st)
        if is_cls
            eps = 1.0f-7
            yc = clamp.(y_hat_val, eps, 1.0f0 - eps)
            per_sample_val = -(y_val .* log.(yc) .+ (1.0f0 .- y_val) .* log.(1.0f0 .- yc))
            current_val = sum(per_sample_val .* reshape(sw_val, 1, :)) / sum(sw_val)
        else
            per_sample_val = (y_hat_val .- y_val) .^ 2
            current_val = sum(per_sample_val .* reshape(sw_val, 1, :)) / sum(sw_val)
        end
        push!(val_losses, current_val)

        if verbose && (epoch % 10 == 0 || epoch == 1)
            @info "Epoch $epoch/$epochs  " *
                  "train=$(round(train_losses[end]; digits=4))  " *
                  "val=$(round(val_losses[end]; digits=4))"
        end

        # Early stopping check
        if use_early_stop
            if current_val < best_val_loss - Float32(min_delta)
                best_val_loss = current_val
                best_ps = deepcopy(ps)
                wait = 0
            else
                wait += 1
                if wait >= patience
                    if verbose
                        @info "Early stopping at epoch $epoch (patience=$patience)"
                    end
                    ps = best_ps
                    break
                end
            end
        end
    end

    model_key = "vnn_$(hash(time()))"
    MODEL_STORE[model_key] = Dict("model" => model, "ps" => ps, "st" => st)

    return Dict(
        "train_losses" => collect(train_losses),
        "val_losses"   => collect(val_losses),
        "model_ref"    => model_key,
    )
end
