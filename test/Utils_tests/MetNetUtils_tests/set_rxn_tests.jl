function test_set_rxn()
    model = deserialize(METNET_CACHE_FILE)

    rxn = Chemostat.Utils.Rxn("new_rxn"; 
        lb = -rand(), ub = rand(), c = rand())
    for met in model.mets
        push!(rxn.mets, met)
        push!(rxn.S, rand([0,-1, 1]))
    end

    rxni = 1
    Chemostat.Utils.set_rxn!(model, rxni, rxn)

    # Test
    @test model.rxns[rxni] == rxn.id
    @test model.lb[rxni] == rxn.lb
    @test model.ub[rxni] == rxn.ub
    @test model.c[rxni] == rxn.c

    for (met, s) in zip(rxn.mets, rxn.S)
        meti = Chemostat.Utils.metindex(model, met)
        @test model.S[meti, rxni] == s
    end
end
test_set_rxn()