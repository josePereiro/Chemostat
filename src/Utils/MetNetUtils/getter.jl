mets(metnet::MetNet) = metnet.mets
mets(metnet::MetNet, ider) = metnet.mets[metindex(metnet, ider)]

metscount(model) = length(model.mets)

rxns(metnet::MetNet) = metnet.rxns
rxns(metnet::MetNet, ider) = metnet.rxns[rxnindex(metnet, ider)]

rxnscount(model) = length(model.rxns)

ub(metnet::MetNet) = metnet.ub
ub(metnet::MetNet, ider) = metnet.ub[rxnindex(metnet, ider)]

lb(metnet::MetNet) = metnet.lb
lb(metnet::MetNet, ider) = metnet.lb[rxnindex(metnet, ider)]

bounds(metnet::MetNet, ider) = (idx = rxnindex(metnet, ider); (metnet.lb[idx], metnet.ub[idx]))

S(metnet::MetNet, metider, rxnider) = 
    metnet.S[metindex(metnet, metider), rxnindex(metnet, rxnider)]

rxn_mets(metnet::MetNet, ider) = findall(metnet.S[:,rxnindex(metnet, ider)] .!= 0.0)
rxn_reacts(metnet::MetNet, ider) = findall(metnet.S[:,rxnindex(metnet, ider)] .< 0.0)
rxn_prods(metnet::MetNet, ider) = findall(metnet.S[:,rxnindex(metnet, ider)] .> 0.0)

met_rxns(metnet::MetNet, ider) = findall(metnet.S[metindex(metnet, ider), :] .!= 0.0)