const fwd_prefix = "_fwd"
const bkwd_prefix = "_bkwd"

using SparseArrays
"""
    Returns a new MetNet with no reversible reactions.
"""
function split_revs(metnet::MetNet, verbose = false) # TODO Add tests
    M, N = size(metnet)
    revs = [isrev(metnet, i) for i in 1:N]
    
    M_, N_ = M, N + count(revs)
    S_ = spzeros(M_, N_)
    S_[1:M,1:N] .= metnet.S
    c_ = copy(metnet.c)
    lb_ = copy(metnet.lb)
    ub_ = copy(metnet.ub)
    grRules_ = copy(metnet.grRules)
    rxns_ = copy(metnet.rxns)
    rxnNames_ = copy(metnet.rxnNames)
    rev_ = falses(N_)
    subSystems_ = copy(metnet.subSystems)
    
    # I only need to update the rev reactions
    _check_and_push!(v::AbstractVector, idx, val) = checkbounds(Bool, v, idx) && push!(v, val)

    for (i, fwd_idx) in enumerate(findall(revs))
        bkwd_idx = N + i
        
        # backup
        orig_rxn = metnet.rxns[fwd_idx]
        orig_lb = metnet.lb[fwd_idx]
        
        # transform to foward reaction
        rxns_[fwd_idx] = "$(orig_rxn)$(fwd_prefix)"
        lb_[fwd_idx] = 0.0
        
        # add backward reaction
        S_[:,bkwd_idx] .= -S_[:,fwd_idx]
        push!(rxns_, "$(orig_rxn)$(bkwd_prefix)")
        push!(rxnNames_, rxnNames_[fwd_idx])
        push!(lb_, 0.0)
        push!(ub_, -orig_lb)
        
        _check_and_push!(c_, bkwd_idx - 1, 0.0) # The objective, if splitted, will be the fwd reaction,
        _check_and_push!(grRules_, bkwd_idx - 1, grRules_[fwd_idx])
        _check_and_push!(subSystems_, bkwd_idx - 1, subSystems_[fwd_idx])

    end

    verbose && N_ != N && (@warn("Only use this method after MetabolicEP.preprocess!!!"); flush(stdout))
    
    return MetNet(metnet; S = S_,c = c_,lb = lb_, ub = ub_, grRules = grRules_, rxns = rxns_, 
        rxnNames =  rxnNames_, rev = rev_, subSystems = subSystems_, 
        intake_info = Dict())
end