function invert_rxn!(metnet::MetNet, ider::IDER_TYPE;
        rename::Union{String,Function} = (rxn) -> string(rxn)
    )
    idx = rxnindex(metnet, ider)
    metnet.S[:,idx] .*= -1
    ub_ = -metnet.lb[idx]
    lb_ = -metnet.ub[idx]
    metnet.ub[idx] = ub_
    metnet.lb[idx] = lb_
    metnet.rxns[idx] = rename isa Function ? rename(metnet.rxns[idx]) : rename
    return metnet
end