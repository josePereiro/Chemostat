function fix_dims(metnet)
    
    function _similar_copy(col, fill, expected_len)
        curr_len = length(col)
        if curr_len < expected_len
            scol = similar(col, expected_len)
            scol[curr_len:end] .= fill
            scol[begin:curr_len] = col
        elseif curr_len > expected_len
            scol = similar(col, expected_len)
            scol[begin:expected_len] = col[begin:expected_len]
        else
            scol = copy(col)
        end
        return scol
    end

    M, N = size(metnet)
    net = Dict()
    net[:S] = copy(metnet.S)
    net[:b] = _similar_copy(metnet.b, 0.0, M)
    net[:c] = _similar_copy(metnet.c, 0.0, N)
    net[:lb] = _similar_copy(metnet.lb, 0.0, N)
    net[:ub] = _similar_copy(metnet.ub, 0.0, N)
    net[:rxns] = copy(metnet.rxns)
    net[:mets] = copy(metnet.mets)
    
    for (l, k) in [(N, :rxns), (M, :mets)]
        length(net[k]) != l && error("Unsolvable miss-match $k ($(length(net[k]))) != $l")
    end

    return MetNet(metnet; reshape = false, net...)
end