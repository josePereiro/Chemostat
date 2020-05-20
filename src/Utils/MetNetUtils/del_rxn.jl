function del_rxn(model::MetNet, iders::Vector)
    to_del_idxs = [rxnindex(model, ider) for ider in iders]
    M, N = size(model)
    left_idxs = filter(x -> !(x in to_del_idxs), 1:N)
    
    _check_and_get(v::AbstractVector) = length(v) == N ? view(v, left_idxs) : v
    
    return MetNet(model; 
                S = view(model.S, 1:M, left_idxs), c = _check_and_get(model.c),
                lb = view(model.lb, left_idxs), ub = view(model.ub, left_idxs),
                genes = _check_and_get(model.genes), grRules = _check_and_get(model.grRules),
                rxns = view(model.rxns, left_idxs), metNames = _check_and_get(model.metNames),
                rev = _check_and_get(model.rev), subSystems = _check_and_get(model.subSystems))
end

del_rxn(model::MetNet, ider::IDER_TYPE) = del_rxn(model::MetNet, [ider])