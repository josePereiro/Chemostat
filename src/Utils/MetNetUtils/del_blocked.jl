function del_blocked(S, b, lb, ub, rxns; eps = 0.0, protected = [])

    lb, ub = (lb, ub) .|> copy # local copy
    m, n = size(S)

    _bidx = trues(n)
    _bidx[protected] .= false
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

function del_blocked(model::MetNet; eps = 0.0, protected = [])
    S, b, lb, ub, rxns = model.S, model.b, model.lb, model.ub, model.rxns
    return del_blocked(S, b, lb, ub, rxns; eps = eps, protected = protected)
end