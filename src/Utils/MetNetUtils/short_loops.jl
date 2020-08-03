function short_loops(model::MetNet; 
        verbose = true)
        
    loops_ = []
    similars_ = similar_rxns(model, verbose = false)

    for(r1, r2) in similars_
        isblock(model, r1) && continue
        isblock(model, r2) && continue
        isfwd(model, r1) && isbkwd(model, r2) && continue
        isbkwd(model, r1) && isfwd(model, r2) && continue

        push!(loops_, (r1, r2))
                
        if verbose
            println(model.rxns[r1], " === ", model.rxns[r2])
            println(model.rxns[r1],": ",rxn_str(model, r1))
            println(model.rxns[r2],": ",rxn_str(model, r2), "\n")
        end
    end

    return loops_
end