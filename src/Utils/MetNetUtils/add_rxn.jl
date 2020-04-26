#TODO fully implement this!!!
function add_rxn(metnet::MetNet, rxnid::AbstractString; kwargs...)
    rxnid in metnet.rxns && error("rxn ($rxnid) already in the model!!!")
    kwargs = Dict(kwargs)
    
    M, N = size(metnet)
    M_, N_ = M, N + 1
    S_ = spzeros(M_, N_)
    S_[1:M,1:N] .= metnet.S
    b_ = copy(metnet.b)
    c_ = copy(metnet.c)
    lb_ = copy(metnet.lb)
    ub_ = copy(metnet.ub)
    genes_ = copy(metnet.genes)
    rxnGeneMat_ = copy(metnet.rxnGeneMat)
    grRules_ = copy(metnet.grRules)
    mets_ = copy(metnet.mets)
    rxns_ = copy(metnet.rxns)
    metNames_ = copy(metnet.metNames)
    metFormulas_ = copy(metnet.metFormulas)
    rxnNames_ = copy(metnet.rxnNames)
    rev_ = copy(metnet.rev)
    subSystems_ = copy(metnet.subSystems)
    
    mets_info = Dict(get(kwargs, :mets, Dict()))
    for (met, stoi_coe) in mets_info
        metidx = metindex(metnet, met)
        S_[metidx, end] = stoi_coe
    end
    
    push!(rxns_, rxnid)
    push!(rxnNames_, get(kwargs, :rxnName, rxnNames_default))
    push!(c_, get(kwargs, :c, c_default))
    push!(lb_, get(kwargs, :lb, lb_default))
    push!(ub_, get(kwargs, :ub, ub_default))
    push!(grRules_, get(kwargs, :grRule, grRules_default))
    push!(subSystems_, get(kwargs, :subSystem, subSystems_default))
    push!(rev_, get(kwargs, :rev, (lb_[end] < 0.0) && (ub_[end] > 0.0)))
    
    return MetNet(N_,M_,S_,b_, c_, lb_, ub_, 
        genes_, rxnGeneMat_, grRules_, mets_, rxns_, 
        metNames_, metFormulas_, rxnNames_, rev_, subSystems_)
    
end