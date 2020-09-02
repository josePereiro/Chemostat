const EMPTY_SPOT = ""
# This just prepare the metnet to hold more elements by making 
# all the numerical fields larger
function expanded_model(metnet::MetNet, newM::Int, newN::Int)
    M, N = size(metnet)
    @assert all((newM, newN) .> (M, N))

    function setfirsts!(col1, col2)
        L = min(length(col1), length(col2))
        col1[1:L] .= col2[1:L]
        return col1
    end
    
    net = Dict()
    net[:S] = zeros(newM, newN)
    net[:S][1:M, 1:N] .= metnet.S
    net[:b] = setfirsts!(zeros(newM), metnet.b)
    net[:c] = setfirsts!(zeros(newN), metnet.c)
    net[:lb] = setfirsts!(zeros(newN), metnet.lb)
    net[:ub] = setfirsts!(zeros(newN), metnet.ub)
    net[:rxns] = setfirsts!(fill(EMPTY_SPOT, newN), metnet.rxns)
    net[:mets] = setfirsts!(fill(EMPTY_SPOT, newM), metnet.mets)
    
    return MetNet(metnet; reshape = false, net...)
end

function findempty(metnet::MetNet, field::Symbol)
    spot = findfirst(isequal(EMPTY_SPOT), getfield(metnet, field))
    isnothing(spot) && error("Not $field empty spot found!!!")
    return spot
end

function compacted_model(metnet::MetNet)

    empty_mets = findall(metnet.mets .== EMPTY_SPOT)
    met_idxs = trues(size(metnet, 1))
    met_idxs[empty_mets] .= false
    M = count(met_idxs)

    empty_rxns = findall(metnet.rxns .== EMPTY_SPOT)
    rxn_idxs = trues(size(metnet, 2))
    rxn_idxs[empty_rxns] .= false
    N = count(rxn_idxs)
    
    net = Dict()
    net[:S] = metnet.S[met_idxs, rxn_idxs]
    net[:b] = metnet.b[met_idxs]
    net[:c] = metnet.c[rxn_idxs]
    net[:lb] = metnet.lb[rxn_idxs]
    net[:ub] = metnet.ub[rxn_idxs]
    net[:rxns] = metnet.rxns[rxn_idxs]
    net[:mets] = metnet.mets[met_idxs]
    
    return MetNet(metnet; net...)
    
end