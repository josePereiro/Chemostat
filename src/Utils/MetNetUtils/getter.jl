mets(metnet::MetNet) = metnet.mets
mets(metnet::MetNet, ider) = metnet.mets[metindex(metnet, ider)]

metscount(model) = length(model.mets)

rxns(metnet::MetNet) = metnet.rxns
rxns(metnet::MetNet, ider) = metnet.rxns[rxnindex(metnet, ider)]

rxnscount(model) = length(model.rxns)


ub(metnet::MetNet, ider) = metnet.ub[rxnindex(metnet, ider)]
lb(metnet::MetNet, ider) = metnet.lb[rxnindex(metnet, ider)]