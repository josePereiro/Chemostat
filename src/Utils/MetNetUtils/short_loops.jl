function short_loops(model::MetNet; 
        verbose = true)
    loops_ = []
    for r1 in eachindex(model.rxns)
        isblock(model, r1) && continue

        cr_reacts = Set(rxn_reacts(model, r1))
        cr_prods = Set(rxn_prods(model, r1))

        for r2 in eachindex(model.rxns)
            isblock(model, r2) && continue
            isfwd(model, r1) && isbkwd(model, r2) && continue
            isbkwd(model, r1) && isfwd(model, r2) && continue

            r_reacts = Set(rxn_reacts(model, r2))
            r_prods = Set(rxn_prods(model, r2))
            if (cr_reacts == r_prods) && (r_reacts == cr_prods) 
                push!(loops_, (r1, r2))
                
                if verbose
                    println(model.rxns[r1], " === ", model.rxns[r2])
                    println(model.rxns[r1],": ",rxn_str(model, r1))
                    println(model.rxns[r2],": ",rxn_str(model, r2), "\n")
                end
            end
        end
    end
    return loops_
end