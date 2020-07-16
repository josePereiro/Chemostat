function del_rxn(model::MetNet, iders::Vector)
    to_del_idxs = [rxnindex(model, ider) for ider in iders]
    M, N = size(model)
    left_idxs = filter(x -> !(x in to_del_idxs), 1:N)
    
    _check_and_get(v::AbstractVector) = length(v) == N ? view(v, left_idxs) : v

    # intake info
    intake_info_ = deepcopy(model.intake_info)
    foreach(rxn -> delete!(intake_info_, rxn), model.rxns[to_del_idxs])
    
    return MetNet(model; 
                S = view(model.S, 1:M, left_idxs), c = _check_and_get(model.c),
                lb = view(model.lb, left_idxs), ub = view(model.ub, left_idxs),
                genes = _check_and_get(model.genes), grRules = _check_and_get(model.grRules),
                rxns = view(model.rxns, left_idxs), rxnNames = _check_and_get(model.rxnNames),
                rev = _check_and_get(model.rev), subSystems = _check_and_get(model.subSystems), 
                intake_info = intake_info_)
end

del_rxn(model::MetNet, ider::IDER_TYPE) = del_rxn(model::MetNet, [ider])