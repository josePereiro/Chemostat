function similar_rxns(model::MetNet; 
    verbose = true)

_cache = Dict()
M, N = size(model)

# find similar reactions
similars_ = []
for r1 in 1:N - 1
    
    if haskey(_cache, r1)
        cr_reacts, cr_prods = _cache[r1]
    else
        cr_reacts = Set(rxn_reacts(model, r1))
        cr_prods = Set(rxn_prods(model, r1))
        _cache[r1] = (cr_reacts, cr_prods)
    end

    for r2 in r1 + 1:N

        if haskey(_cache, r2)
            r_reacts, r_prods = _cache[r2]
        else
            r_reacts = Set(rxn_reacts(model, r2))
            r_prods = Set(rxn_prods(model, r2))
            _cache[r2] = (r_reacts, r_prods)
        end
        
        if (cr_reacts == r_prods) && (r_reacts == cr_prods) 
            push!(similars_, (r1, r2))
            
            if verbose
                println(model.rxns[r1], " === ", model.rxns[r2])
                println(model.rxns[r1], ": ", rxn_str(model, r1))
                println(model.rxns[r2], ": ", rxn_str(model, r2), "\n")
            end
        end
    end
end
return similars_
end