#TODO fully implement this!!!
# TODO: vectorize this
# TODO: make vectorized version
"""
    Will add a Met to a given idx. See Met struct for details
"""
function set_rxn!(metnet::MetNet, rxni::Int, rxn::Rxn)
    metnet.rxns[rxni] = rxn.id
    metnet.lb[rxni] = rxn.lb
    metnet.ub[rxni] = rxn.ub
    metnet.c[rxni] = rxn.c
    for (met, s) in zip(rxn.mets, rxn.S)
        meti = metindex(metnet, met)
        metnet.S[meti, rxni] = s
    end
    return metnet
end

set_rxn!(metnet::MetNet, rxn::Rxn) = 
    set_rxn!(metnet, findempty(metnet, :rxns; check = true), rxn)
    
set_rxn!(metnet::MetNet, rxni::Int, id::String) = set_rxn!(metnet, rxni, Rxn(id))
set_rxn!(metnet::MetNet, id::String) = 
    set_rxn!(metnet, findempty(metnet, :rxns; check = true), id)