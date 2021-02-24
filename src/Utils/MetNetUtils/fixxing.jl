function fixxing(f::Function, model::MetNet, ider, val::Real)
    idx = rxnindex(model, ider)
    bk_lb, bk_ub = model.lb[idx], model.ub[idx]
    model.lb[idx] = model.ub[idx] = val
    try
        val = f()
    finally
        model.lb[idx], model.ub[idx] = bk_lb, bk_ub
    end
    return val
end
