"""
    Findall the reaction that transforms the same elements
    if we define: r1 = 2A -> 2B, r2 = B <-> A, r3 = A <- B, r4 = A + B -> 
    the method will returns [[r1, r2, r3]]
    This method only use the information in the binary S. 
    Bounds are not relevant either.
"""
function similar_rxns(model::MetNet, 
        iders = eachindex(model.rxns); verbose = true)

    idxs = [rxnindex(model, ider) for ider in iders]
    # collecting react and prods hashs
    # (prods_hash, react_hash) => reaction idxs
    # (react_hash, prods_hash) => reaction idxs
    hash_table = Dict{Tuple{UInt64,UInt64}, Vector{Int}}()
    for rxni in idxs
        rhash = hash(Set(rxn_reacts(model, rxni)))
        phash = hash(Set(rxn_prods(model, rxni)))
        rxnis = get!(hash_table, (rhash, phash), Int[])
        push!(rxnis, rxni)
    end

    # Form pairs
    similars = Vector{Vector{Int}}()
    for (fkey, frxnis) in hash_table

        bkey = fkey |> reverse
        brxnis = get(hash_table, bkey, Int[])

        sims_ = union(frxnis, brxnis)
        length(sims_) > 1 && push!(similars, sims_)

        # deleting
        delete!(hash_table, fkey)
        delete!(hash_table, bkey)

    end

    if verbose
        for rxnis in similars
            println("Similars ", length(rxnis))
            for rxni in rxnis
                println(model.rxns[rxni], " [", rxni, "]: ", rxn_str(model, rxni))
            end
            println()
        end
    end
    return similars
end

