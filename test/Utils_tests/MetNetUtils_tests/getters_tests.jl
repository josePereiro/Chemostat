function test_getters()
    model = deserialize(METNET_CACHE_FILE)
    @test Chemostat.Utils.mets(model) == model.mets
    @test Chemostat.Utils.rxns(model) == model.rxns

    for (i, rxn) in enumerate(Chemostat.Utils.rxns(model))
        @test Chemostat.Utils.ub(model, rxn) == model.ub[i]
        @test Chemostat.Utils.ub(model, i) == model.ub[i]
        @test Chemostat.Utils.lb(model, rxn) == model.lb[i]
        @test Chemostat.Utils.lb(model, i) == model.lb[i]
        @test Chemostat.Utils.rxns(model, i) == model.rxns[i]
        @test Chemostat.Utils.rxns(model, rxn) == rxn
    end

    for (i, met) in enumerate(Chemostat.Utils.mets(model))
        @test Chemostat.Utils.mets(model, i) == model.mets[i]
        @test Chemostat.Utils.mets(model, met) == model.mets[i]
    end

end
test_getters()