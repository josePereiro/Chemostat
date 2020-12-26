"""
    Will add a Met to a given idx. See Met struct for details
"""
function set_met!(metnet::MetNet, meti::Int, met::Met)
    metnet.mets[meti] = met.id
    for (rxn, s) in zip(met.rxns, met.S)
        rxni = rxnindex(metnet, rxn)
        metnet.S[meti, rxni] = s
    end
    metnet.b[meti] = met.b
    return metnet
end

set_met!(metnet::MetNet, met::Met) = 
    set_met!(metnet, findempty(metnet, :mets; check = true), met)

set_met!(metnet::MetNet, meti::Int, id::String) = set_met!(metnet, meti, Met(id))
set_met!(metnet::MetNet, id::String) = 
    set_met!(metnet, findempty(metnet, :mets; check = true), id)