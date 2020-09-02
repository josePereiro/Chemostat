function del_met!(metnet::MetNet, met::IDER_TYPE)
    meti = metindex(metnet, met)
    metnet.mets[meti] = EMPTY_SPOT
    z = zero(eltype(metnet.S))
    metnet.S[meti, :] .= z
    metnet.b[meti] = z
    
    return metnet
end