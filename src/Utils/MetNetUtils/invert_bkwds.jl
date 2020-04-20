function invert_bkwds!(metnet)
    revscount(metnet) > 0 && error("All reactions must be irreversibles!!!")

    for bkwd_rxn in bkwds(metnet)
        metnet.S[:,bkwd_rxn] .*= -1
        lb_ = abs(metnet.ub[bkwd_rxn])
        ub_ = abs(metnet.lb[bkwd_rxn])
        metnet.ub[bkwd_rxn] = ub_
        metnet.lb[bkwd_rxn] = lb_
        metnet.rxns[bkwd_rxn] *= bkwd_prefix
    end
    
    return metnet
end
invert_bkwds(metnet) = invert_bkwds!(deepcopy(metnet))