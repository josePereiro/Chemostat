function del_met(model::MetNet, iders::Vector)

    to_del_idxs = [metindex(model, ider) for ider in iders]
    M, N = size(model)

    bidx = trues(M)
    bidx[to_del_idxs] .= false
    left_idxs = collect(1:M)[bidx]
    
    _check_and_get(v::AbstractVector) = length(v) == M ? view(v, left_idxs) : v
    
    return MetNet(model; 
                S = view(model.S, left_idxs, 1:N),
                b = view(model.b, left_idxs),
                mets = view(model.mets, left_idxs), 
                metNames = _check_and_get(model.metNames), 
                metFormulas = _check_and_get(model.metFormulas))
end

del_met(model::MetNet, ider::IDER_TYPE) = del_met(model::MetNet, [ider])