function del_blocked(S, b, lb, ub, rxns; eps = 0.0, protect = [])

    lb, ub = (lb, ub) .|> copy # local copy
    m, n = size(S)

    _bidx = trues(n)
    _bidx[protect] .= false
    non_protected = findall(_bidx)

    blocked = falses(n)
    for i in non_protected
        if lb[i] == ub[i] # blocked
            blocked[i] = eps == 0.0
            lb[i], ub[i] = lb[i] - eps, ub[i] + eps
        end
    end 
    
    unblocked = findall(.!blocked)

    return S[:,unblocked], b - S*(blocked .* lb), 
        lb[unblocked], ub[unblocked], rxns[unblocked], findall(blocked)
end

function del_blocked(model::MetNet; eps = 0.0, protect = [])
    protect = map((r) -> rxnindex(model, r), protect)
    S, b, lb, ub, rxns, blocked = 
        del_blocked(model.S, model.b, model.lb, model.ub, model.rxns; eps, protect)
    return MetNet(model; reshape = true, S, b, lb, ub, rxns)
end