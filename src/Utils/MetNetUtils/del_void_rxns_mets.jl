function del_void_rxns_mets!(metnet::MetNet)
    # find mets with no rxns
    for (meti, met) in enumerate(metnet.mets)
        rxns = met_rxns(metnet, meti)
        isempty(rxns) && del_met!(metnet, meti)
    end

    # find rxns with no mets
    for (rxni, rxn) in enumerate(metnet.rxns)
        mets = rxn_mets(metnet, rxni)
        isempty(mets) && del_rxn!(metnet, rxni)
    end

    return metnet
end