mets(metnet::MetNet) = metnet.mets
mets(metnet::MetNet, ider::IDER_TYPE) = metnet.mets[metindex(metnet, ider)]

metscount(metnet::MetNet) = length(metnet.mets)

rxns(metnet::MetNet) = metnet.rxns
rxns(metnet::MetNet, ider::IDER_TYPE) = metnet.rxns[rxnindex(metnet, ider)]

rxnscount(model) = length(model.rxns)

ub(metnet::MetNet) = metnet.ub
ub(metnet::MetNet, ider::IDER_TYPE) = metnet.ub[rxnindex(metnet, ider)]

lb(metnet::MetNet) = metnet.lb
lb(metnet::MetNet, ider::IDER_TYPE) = metnet.lb[rxnindex(metnet, ider)]

b(metnet::MetNet) = metnet.b
b(metnet::MetNet, ider::IDER_TYPE) = metnet.b[rxnindex(metnet, ider)]

bounds(metnet::MetNet, ider::IDER_TYPE) = (idx = rxnindex(metnet, ider); (metnet.lb[idx], metnet.ub[idx]))

S(metnet::MetNet, metider::IDER_TYPE, rxnider::IDER_TYPE) = 
    metnet.S[metindex(metnet, metider), rxnindex(metnet, rxnider)]

rxn_mets(metnet::MetNet, ider::IDER_TYPE) = findall(metnet.S[:,rxnindex(metnet, ider)] .!= 0.0)
rxn_reacts(metnet::MetNet, ider::IDER_TYPE) = findall(metnet.S[:,rxnindex(metnet, ider)] .< 0.0)
rxn_prods(metnet::MetNet, ider::IDER_TYPE) = findall(metnet.S[:,rxnindex(metnet, ider)] .> 0.0)

met_rxns(metnet::MetNet, ider::IDER_TYPE) = findall(metnet.S[metindex(metnet, ider), :] .!= 0.0)