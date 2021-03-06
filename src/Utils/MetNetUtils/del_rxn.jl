function del_rxn!(metnet::MetNet, rxn::IDER_TYPE)
    rxni = rxnindex(metnet, rxn)
    metnet.rxns[rxni] = EMPTY_SPOT
    z = zero(eltype(metnet.S))
    metnet.S[:, rxni] .= z
    metnet.lb[rxni] = z
    metnet.ub[rxni] = z
    metnet.c[rxni] = z
    
    return metnet
end

function del_rxn!(metnet::MetNet, rxns::Vector)
    for rxn in rxns
        del_rxn!(metnet, rxn)
    end
    return metnet
end