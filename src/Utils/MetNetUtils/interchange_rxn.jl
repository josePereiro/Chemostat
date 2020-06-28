# TODO fully implement this function
function interchange_rxn!(model::MetNet, rxn1_ider::IDER_TYPE, rxn2_ider::IDER_TYPE)
    rxn1_idx = rxnindex(model, rxn1_ider)
    rxn2_idx = rxnindex(model, rxn2_ider)
    
    # S
    temp = model.S[:, rxn1_idx]
    model.S[:, rxn1_idx] .= model.S[:, rxn2_idx]
    model.S[:, rxn2_idx] .= temp
    
    # rxn
    temp = model.rxns[rxn1_idx]
    model.rxns[rxn1_idx] = model.rxns[rxn2_idx]
    model.rxns[rxn2_idx] = temp
    
    # rxnNames
    temp = model.rxnNames[rxn1_idx]
    model.rxnNames[rxn1_idx] = model.rxnNames[rxn2_idx]
    model.rxnNames[rxn2_idx] = temp
    
    # lb
    temp = model.lb[rxn1_idx]
    model.lb[rxn1_idx] = model.lb[rxn2_idx]
    model.lb[rxn2_idx] = temp
    
    # ub
    temp = model.ub[rxn1_idx]
    model.ub[rxn1_idx] = model.ub[rxn2_idx]
    model.ub[rxn2_idx] = temp
    
    # rev
    temp = model.rev[rxn1_idx]
    model.rev[rxn1_idx] = model.rev[rxn2_idx]
    model.rev[rxn2_idx] = temp
    
end