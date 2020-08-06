"""
    Returns a metnet with the added metabolite metid.

    kwargs
    - rxns: a key-pair data type with the ider of the reactions and the stoichiometric coefficient
        associated with metid
    Ex: rxns = Dict("rxn1" => -1, "rxn2" => 1)
    - metName: the long name of the metabolite
    - metFormula: what do you think?
"""
function add_met(metnet::MetNet, metid::AbstractString; kwargs...)
    metid in metnet.mets && error("met ($metid) already in the model!!!")
    kwargs = Dict(kwargs)
    
    M, N = size(metnet)
    S_ = spzeros(M + 1, N)
    S_[1:M,1:N] .= metnet.S
    rxns_info = Dict(get(kwargs, :rxns, Dict()))
    for (rxn, stoi_coe) in rxns_info
        rxnidx = rxnindex(metnet, rxn)
        S_[end, rxnidx] = stoi_coe
    end

    b_ = vcat(metnet.b, get(kwargs, :b, b_default))
    mets_ = vcat(metnet.mets, metid)
    metNames_ = vcat(metnet.metNames, get(kwargs, :metName, metNames_default))
    metFormulas_ = vcat(metnet.metFormulas, get(kwargs, :metFormula, metFormulas_default))
    
    return MetNet(metnet; reshape = true, S = S_, b = b_, mets = mets_, 
                    metNames = metNames_, metFormulas = metFormulas_)
end