function clamp_bounds!(model, max_abs = 1e3, zeroth = 1e-8)
    M, N = size(model)
    @inbounds for i in 1:N
        lb = model.lb[i]
        ub = model.ub[i]
        model.lb[i] = clamp(abs(lb) < zeroth ? zero(lb) : lb, -max_abs, max_abs)
        model.ub[i] = clamp(abs(ub) < zeroth ? zero(ub) : ub, -max_abs, max_abs)
    end
    return model
end