function invert_bkwds!(metnet::MetNet; rename = true)
    for bkwd_rxn in bkwds_bounded(metnet)
        invert_rxn!(metnet, bkwd_rxn, rename = rename)
    end
    return metnet
end
invert_bkwds(metnet::MetNet) = invert_bkwds!(deepcopy(metnet))