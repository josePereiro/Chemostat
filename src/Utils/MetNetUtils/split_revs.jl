const fwd_prefix = "_fwd"
const bkwd_prefix = "_bkwd"

using SparseArrays
"""
    Returns a new MetNet with no reversible reactions.
"""
function split_revs(metnet::MetNet) # TODO Add tests
    M, N = metnet.M, metnet.N
    revs = [isrev(metnet, i) for i in 1:N]
    
    M_, N_ = M, N + count(revs)
    S_ = spzeros(M_, N_)
    S_[1:M,1:N] .= metnet.S
    c_ = copy(metnet.c)
    lb_ = copy(metnet.lb)
    ub_ = copy(metnet.ub)
    genes_ = copy(metnet.genes)
    grRules_ = copy(metnet.grRules)
    rxns_ = copy(metnet.rxns)
    rxnNames_ = copy(metnet.rxnNames)
    rev_ = falses(N_)
    subSystems_ = copy(metnet.subSystems)
    
    # I only need to update the rev reactions
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

        checkbounds(Bool, c_, bkwd_idx - 1) && push!(c_, 0.0) # The objective, if splitted, will be the fwd reaction,
        checkbounds(Bool, grRules_, bkwd_idx - 1) && push!(grRules_, grRules_[fwd_idx])
        checkbounds(Bool, subSystems_, bkwd_idx - 1) && push!(subSystems_, subSystems_[fwd_idx])
        checkbounds(Bool, genes_, bkwd_idx - 1) && push!(genes_, genes_[fwd_idx])

    end

    N_ != N && (@warn("Only use this method after MetabolicEP.preprocess!!!"); flush(stdout))
    
    return MetNet(S = S_,c = c_,lb = lb_, ub = ub_, 
        genes = genes_, grRules = grRules_, rxns = rxns_, 
        rxnNames =  rxnNames_, rev = rev_, subSystems = subSystems_)
end