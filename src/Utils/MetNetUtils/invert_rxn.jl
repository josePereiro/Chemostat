function invert_rxn!(metnet::MetNet, ider::IDER_TYPE;
        rename = true)
    idx = rxnindex(metnet, ider)
    metnet.S[:,idx] .*= -1
    lb_ = abs(metnet.ub[idx])
    ub_ = abs(metnet.lb[idx])
    metnet.ub[idx] = ub_
    metnet.lb[idx] = lb_
    rename && (metnet.rxns[idx] *= bkwd_prefix) 
    return metnet
end