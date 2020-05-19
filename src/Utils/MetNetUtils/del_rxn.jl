function del_rxn(model, ider)
    rxn_idx = rxnindex(model, ider)
    N = size(model, 2)
    left_idxs = 1:N .!= rxn_idx
    
    _check_and_get(v::AbstractVector) = length(v) == N ? view(v, left_idxs) : v
    
    return MetNet(model; S = view(model.S, :, left_idxs), 
        c = _check_and_get(model.c),
        lb = view(model.lb, left_idxs),
        ub = view(model.ub, left_idxs),
        genes = _check_and_get(model.genes),
        grRules = _check_and_get(model.grRules),
        rxns = view(model.rxns, left_idxs),
        metNames = _check_and_get(model.metNames),
        rev = _check_and_get(model.rev),
        subSystems = _check_and_get(model.subSystems)
    )
end