const REV_SUFIX = "_REV"
function invert_bkwds!(metnet::MetNet; rename::Function = (rxn) -> string(rxn, REV_SUFIX))
    for bkwd_rxn in bkwds_bounded(metnet)
        invert_rxn!(metnet, bkwd_rxn; rename)
    end
    return metnet
end
invert_bkwds(metnet::MetNet; rename::Function = (rxn) -> rxn) = 
    invert_bkwds!(deepcopy(metnet); rename)