function fixxing(f::Function, model::MetNet, ider, val::Real; 
        btol = 0.0
    )
    idx = rxnindex(model, ider)
    bk_lb, bk_ub = model.lb[idx], model.ub[idx]
    try
        dflx = abs(val * btol)
        model.lb[idx] = val - dflx
        model.ub[idx] = val + dflx
        return f()
    finally
        model.lb[idx], model.ub[idx] = bk_lb, bk_ub
    end
end
