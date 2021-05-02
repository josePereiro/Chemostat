const FWD_SUFFIX = "_fwd"
const BKWD_SUFFIX = "_bkwd"

function split_rxns(metnet::MetNet, iders;
        get_fwd_ider::Function = (rxn) -> string(rxn, FWD_SUFFIX),
        get_bkwd_ider::Function = (rxn) -> string(rxn, BKWD_SUFFIX),
        on_split::Function = (new_model, fwd_idx, bkwd_idx) -> nothing
    )

    M, N = size(metnet)
    idxs = rxnindex.([metnet], iders)
    newM, newN = M, N + length(idxs)

    new_model = expanded_model(metnet, newM, newN)

    for (i, fwd_idx) in enumerate(idxs)
        
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
        
        on_split(new_model, fwd_idx, bkwd_idx)
    end

    return new_model

end

"""
    Returns a new MetNet with no reversible reactions.
"""
function split_revs(metnet::MetNet;
        get_fwd_ider::Function = (rxn) -> string(rxn, FWD_SUFFIX),
        get_bkwd_ider::Function = (rxn) -> string(rxn, BKWD_SUFFIX),
        on_rev::Function = (new_model, fwd_idx, bkwd_idx) -> nothing
    ) # TODO Add tests

    return split_rxns(metnet, revs(metnet); 
        get_fwd_ider, get_bkwd_ider, on_split = on_rev
    )
end

# function split_revs(metnet::MetNet;
#         get_fwd_ider::Function = (rxn) -> string(rxn, FWD_SUFFIX),
#         get_bkwd_ider::Function = (rxn) -> string(rxn, BKWD_SUFFIX),
#         on_rev!::Function = (new_model, fwd_idx, bkwd_idx) -> nothing
#     ) # TODO Add tests

#     M, N = size(metnet)
#     revs = revs(metnet)
#     newM, newN = M, N + length(revs)

#     new_model = expanded_model(metnet, newM, newN)
    
#     # I only need to update the rev reactions
#     _check_and_push!(v::AbstractVector, idx, val) = checkbounds(Bool, v, idx) && push!(v, val)

#     for (i, fwd_idx) in enumerate(findall(revs))
        
#         # backup
#         orig_rxn = metnet.rxns[fwd_idx]
#         orig_lb = metnet.lb[fwd_idx]
        
#         # transform to forward reaction
#         new_model.rxns[fwd_idx] = get_fwd_ider(orig_rxn)
#         new_model.lb[fwd_idx] = 0.0
        
#         # add backward reaction
#         bkwd_idx = findempty(new_model, :rxns)
#         new_model.rxns[bkwd_idx] = get_bkwd_ider(orig_rxn)
#         new_model.S[:,bkwd_idx] .= -new_model.S[:,fwd_idx]
#         new_model.c[bkwd_idx] = 0.0
#         new_model.lb[bkwd_idx] = 0.0
#         new_model.ub[bkwd_idx] = abs(orig_lb)
        
#         on_rev!(new_model, fwd_idx, bkwd_idx)
#     end

#     return new_model
# end