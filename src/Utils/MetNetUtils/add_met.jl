# TODO fully implement this!!!
"""
    Returns a metnet with the added metabolite metid.

    kwargs
    - rxns: a key-pair data type with the ider of the reactions and the stoichiometric coefficient
        associated with metid
    Ex: rxns = Dict("rxn1" => -1, "rxn2" => 1)
    - metName: the long name of the metabolite
    - metFormula: what do you think?
"""
function add_met(metnet, metid::AbstractString; kwargs...)
    metid in metnet.mets && error("met ($metid) already in the model!!!")
    kwargs = Dict(kwargs)
    
    M, N = size(metnet)
    M_, N_ = M + 1, N
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
    
    rxns_info = Dict(get(kwargs, :rxns, Dict()))
    for (rxn, stoi_coe) in rxns_info
        rxnidx = rxnindex(metnet, rxn)
        S_[end, rxnidx] = stoi_coe
    end
    
    push!(mets_, metid)
    push!(metNames_, get(kwargs, :metName, ""))
    push!(metFormulas_, get(kwargs, :metFormula, ""))
    
    return MetNet(N_,M_,S_,b_, c_, lb_, ub_, 
        genes_, rxnGeneMat_, grRules_, mets_, rxns_, 
        metNames_, metFormulas_, rxnNames_, rev_, subSystems_)
end