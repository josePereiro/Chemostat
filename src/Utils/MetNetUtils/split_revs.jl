const FWD_SUFFIX = "_fwd"
const BKWD_SUFFIX = "_bkwd"

using SparseArrays
"""
    Returns a new MetNet with no reversible reactions.
"""
function split_revs(metnet::MetNet;
        get_fwd_ider::Function = (rxn) -> string(rxn, FWD_SUFFIX),
        get_bkwd_ider::Function = (rxn) -> string(rxn, BKWD_SUFFIX),
        on_rev!::Function = (new_model, fwd_idx, bkwd_idx) -> nothing
    ) # TODO Add tests

    M, N = size(metnet)
    revs = [isrev(metnet, i) for i in 1:N]
    newM, newN = M, N + count(revs)

    new_model = expanded_model(metnet, newM, newN)
    
    # I only need to update the rev reactions
    _check_and_push!(v::AbstractVector, idx, val) = checkbounds(Bool, v, idx) && push!(v, val)

    for (i, fwd_idx) in enumerate(findall(revs))
        
        # backup
        orig_rxn = metnet.rxns[fwd_idx]
        orig_lb = metnet.lb[fwd_idx]
        
        # transform to forward reaction
        new_model.rxns[fwd_idx] = get_fwd_ider(orig_rxn)
        new_model.lb[fwd_idx] = 0.0
        
        # add backward reaction
        bkwd_idx = findempty(new_model, :rxns)
        new_model.rxns[bkwd_idx] = get_bkwd_ider(orig_rxn)
        new_model.S[:,bkwd_idx] .= -new_model.S[:,fwd_idx]
        new_model.c[bkwd_idx] = 0.0
        new_model.lb[bkwd_idx] = 0.0
        new_model.ub[bkwd_idx] = abs(orig_lb)
        
        on_rev!(new_model, fwd_idx, bkwd_idx)
    end

    return new_model
end