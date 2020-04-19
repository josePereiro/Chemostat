function test_TestUtils_linear_S()
    S =[1.0 -1.0  0.0  0.0;
        0.0  1.0 -1.0  0.0;
        0.0  0.0  1.0 -1.0]
    @test  all(Chemostat.Utils.lineal_S(size(S,1)) .== S)
end
test_TestUtils_linear_S()

function test_MetNetUtils_iders()
    model = Chemostat.Utils.simple_toy_MetNet()

    @test_throws ErrorException Chemostat.Utils.metindex(model, 100)
    @test_throws ErrorException Chemostat.Utils.metindex(model, -100)

    @test_throws ErrorException Chemostat.Utils.rxnindex(model, 100)
    @test_throws ErrorException Chemostat.Utils.rxnindex(model, -100)


    for (i, rxn) in enumerate(Chemostat.Utils.rxns(model))
        @test Chemostat.Utils.rxnindex(model, i) == i
        @test Chemostat.Utils.rxnindex(model, rxn) == i
    end

    for (i, met) in enumerate(Chemostat.Utils.mets(model))
        @test Chemostat.Utils.metindex(model, i) == i
        @test Chemostat.Utils.metindex(model, met) == i
    end

    @test_throws ErrorException Chemostat.Utils.metindex(model, "met_not_in_model")
    @test_throws ErrorException Chemostat.Utils.rxnindex(model, "rxn_not_in_model")

end
test_MetNetUtils_iders()

function test_MetNetUtils_getters()
    model = Chemostat.Utils.simple_toy_MetNet()
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
test_MetNetUtils_getters()

