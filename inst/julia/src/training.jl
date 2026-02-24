# =============================================================================
# training.jl — Training loop with biology-aware regularization + early stopping
# =============================================================================

const MODEL_STORE = Dict{String, Any}()

# Handle single mask name from R
# JuliaConnectoR sends length-1 R character vectors as scalar String
function train_vnn(mask_name::String, args...)
    return train_vnn([mask_name], args...)
end

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
    rng = Random.MersenneTwister(seed)

    X_all = DATA_STORE["X"]
    y_all = DATA_STORE["y"]
    n_output = DATA_STORE["n_output"]
    act_name = DATA_STORE["activation"]

    n_samples = size(X_all, 1)

    # Train/val split
    n_val = max(1, round(Int, n_samples * val_frac))
    perm = randperm(rng, n_samples)
    val_idx = perm[1:n_val]
    train_idx = perm[n_val+1:end]

    X_train = X_all[train_idx, :]'
    y_train = reshape(y_all[train_idx], 1, :)
    X_val   = X_all[val_idx, :]'
    y_val   = reshape(y_all[val_idx], 1, :)
    n_train = length(train_idx)

    # Build model
    act_fn = get_activation_fn(act_name)
    model = build_vnn_model(mask_names; n_output = n_output, activation_fn = act_fn)
    ps, st = Lux.setup(rng, model)

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

            result = Zygote.withgradient(ps) do p
                y_hat, _ = model(xb, p, st)

                if is_cls
                    eps = 1.0f-7
                    yc = clamp.(y_hat, eps, 1.0f0 - eps)
                    data_loss = -mean(yb .* log.(yc) .+ (1.0f0 .- yb) .* log.(1.0f0 .- yc))
                else
                    data_loss = mean((y_hat .- yb) .^ 2)
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

        # Validation loss
        y_hat_val, _ = model(X_val, ps, st)
        if is_cls
            eps = 1.0f-7
            yc = clamp.(y_hat_val, eps, 1.0f0 - eps)
            current_val = -mean(y_val .* log.(yc) .+ (1.0f0 .- y_val) .* log.(1.0f0 .- yc))
        else
            current_val = mean((y_hat_val .- y_val) .^ 2)
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
