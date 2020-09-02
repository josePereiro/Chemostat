function test_set_met()
    model = deserialize(METNET_CACHE_FILE)

    met = Chemostat.Utils.Met("new_met"; b = rand())
    for rxn in model.rxns
        push!(met.rxns, rxn)
        push!(met.S, rand([0,-1, 1]))
    end
    
    meti = 1
    Chemostat.Utils.set_met!(model, meti, met)
    
    # Tests
    @test model.mets[meti] == met.id
    @test model.b[meti] == met.b

    for (rxn, s) in zip(met.rxns, met.S)
        rxni = Chemostat.Utils.rxnindex(model, rxn)
        @test model.S[meti, rxni] == s
    end

end
test_set_met()