function similar_rxns(model::MetNet; 
    verbose = true)

    M, N = size(model)
    # collecting react and prods hashs
    _rxn_hashs = Vector{Tuple{UInt64,UInt64}}(undef, N)
    for r in 1:N
        reacts = hash(Set(rxn_reacts(model, r)))
        prods = hash(Set(rxn_prods(model, r)))
        _rxn_hashs[r] = (reacts, prods)
    end


    # find similar reactions
    similars_ = []
    for r1 in 1:N - 1

        cr_reacts, cr_prods = _rxn_hashs[r1]
        
        for r2 in r1 + 1:N

            r_reacts, r_prods = _rxn_hashs[r2]
            
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