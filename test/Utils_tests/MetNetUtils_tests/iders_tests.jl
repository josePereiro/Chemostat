function test_iders()
    model = Chemostat.Utils.simple_toy_MetNet()

    too_big = 100
    negative = -100
    @test_throws ErrorException Chemostat.Utils.metindex(model, too_big)
    @test_throws ErrorException Chemostat.Utils.metindex(model, negative)

    @test_throws ErrorException Chemostat.Utils.rxnindex(model, too_big)
    @test_throws ErrorException Chemostat.Utils.rxnindex(model, negative)


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
test_iders()