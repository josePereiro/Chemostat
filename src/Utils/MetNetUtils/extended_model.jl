const EMPTY_SPOT = ""
# This just prepare the metnet to hold more elements by making 
# all the numerical fields larger
function expanded_model(metnet::MetNet{T}, newM::Int, newN::Int) where T
    M, N = size(metnet)
    @assert all((newM, newN) .>= (M, N))
    check_dims(metnet)

    function _similar_copy(col, fill, newdim)
        L = length(col)
        newcol = similar(col, newdim)
        newcol[L:end] .= fill
        newcol[1:L] .= col[1:L]
        return newcol
    end

    net = Dict()

    net[:S] = similar(metnet.S, newM, newN)
    net[:S][M:end, :] .= zero(T)
    net[:S][:, N:end] .= zero(T)
    net[:S][1:M, 1:N] .= metnet.S

    net[:b] = _similar_copy(metnet.b, 0.0, newM)
    net[:c] = _similar_copy(metnet.c, 0.0, newN)
    net[:lb] = _similar_copy(metnet.lb, 0.0, newN)
    net[:ub] = _similar_copy(metnet.ub, 0.0, newN)
    net[:rxns] = _similar_copy(metnet.rxns, EMPTY_SPOT, newN)
    net[:mets] = _similar_copy(metnet.mets, EMPTY_SPOT, newM)
    
    return MetNet(metnet; reshape = false, net...)
end

function findempty(metnet::MetNet, field::Symbol; check = false) 
    idx = findfirst(isequal(EMPTY_SPOT), getfield(metnet, field))
    check && isnothing(idx) && error("no empty spot found in $(field) ")
    idx
end

is_compacted(metnet::MetNet) = !isnothing(findempty(metnet, :mets)) && !isnothing(findempty(metnet, :rxns))

function compacted_model(metnet::MetNet)

    empty_mets = findall(metnet.mets .== EMPTY_SPOT)
    met_idxs = trues(size(metnet, 1))
    met_idxs[empty_mets] .= false

    empty_rxns = findall(metnet.rxns .== EMPTY_SPOT)
    rxn_idxs = trues(size(metnet, 2))
    rxn_idxs[empty_rxns] .= false
    
    net = Dict()
    net[:S] = metnet.S[met_idxs, rxn_idxs]
    net[:b] = metnet.b[met_idxs]
    net[:c] = metnet.c[rxn_idxs]
    net[:lb] = metnet.lb[rxn_idxs]
    net[:ub] = metnet.ub[rxn_idxs]
    net[:rxns] = metnet.rxns[rxn_idxs]
    net[:mets] = metnet.mets[met_idxs]
    
    return MetNet(metnet; reshape = false, net...)
    
end