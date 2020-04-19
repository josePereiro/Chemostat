const fwd_prefix = "_fwd"
const bkwd_prefix = "_bfwd"

using SparseArrays
"""
    Retiurns a new MetNet with no reversible reactions.
"""
function split_revs(metnet)
    M,N = metnet.M, metnet.N
    revs = [isrev(metnet, i) for i in 1:N]
    
    M_, N_ = M, N + count(revs)
    S_ = spzeros(M_, N_)
    S_[1:M,1:N] .= metnet.S
    b_ = copy(metnet.b)
    c_ = copy(metnet.c)
    lb_ = copy(metnet.lb)
    ub_ = copy(metnet.ub)
    #TODO
    genes_ = copy(metnet.genes)
    rev_ = Bool.(zeros(N))
    rxnGeneMat_ = copy(metnet.rxnGeneMat)
    grRules_ = copy(metnet.grRules)
    mets_ = copy(metnet.mets)
    rxns_ = copy(metnet.rxns)
    metNames_ = copy(metnet.metNames)
    metFormulas_ = copy(metnet.metFormulas)
    rxnNames_ = copy(metnet.rxnNames)
    revs_ = Bool.(zeros(N_))
    subSystems_ = copy(metnet.subSystems)
    
    # I only need to update the rev reactions
    for (i, fwd_idx) in enumerate(findall(revs))
        bkwd_idx = N + i
        
        # backup
        orig_rxn = metnet.rxns[fwd_idx]
        orig_lb = metnet.lb[fwd_idx]
        
        println(orig_rxn)
        println(i)
        
        
        # foward reaction
        rxns_[fwd_idx] = "$(orig_rxn)$(fwd_prefix)"
        lb_[fwd_idx] = 0.0
        
        # backward reaction
        S_[:,bkwd_idx] .= -S_[:,fwd_idx]
        push!(rxns_, "$(orig_rxn)$(bkwd_prefix)")
        push!(rxnNames_, rxnNames_[fwd_idx])
        push!(grRules_, grRules_[fwd_idx])
        push!(subSystems_, subSystems_[fwd_idx])
        push!(lb_, 0.0)
        push!(ub_, -orig_lb)
        push!(c_, 0.0) # The objective, if splitted, will be the fwd reaction
        
        
    end
    
    return MetNet(N_,M_,S_,b_, c_, lb_, ub_, 
        genes_, rxnGeneMat_, grRules_, mets_, rxns_, 
        metNames_, metFormulas_, rxnNames_, rev_, subSystems_)
end