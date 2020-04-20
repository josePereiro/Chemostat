function test_MetNetUtils_iders()
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

function test_MetNetUtils_add_met(model, rxns_info_)
    
    metid_ = "new_met"
    metName_ = "juanito"
    metFormula_ = "1+1=3"
    new_model = Chemostat.Utils.add_met(model, metid_, 
        rxns = rxns_info_, metName = metName_, metFormula = metFormula_)

    @test new_model.mets[end] == metid_
    @test new_model.metNames[end] == metName_
    @test new_model.metFormulas[end] == metFormula_
    @test_throws ErrorException Chemostat.Utils.add_met(model, model.mets[1])

    for (rxn, stoi_coe) in rxns_info_
        rxn_indx = Chemostat.Utils.rxnindex(new_model, rxn)
        @test new_model.S[end, rxn_indx] == stoi_coe
    end

end

function test_MetNetUtils_add_met()
    model = Chemostat.Utils.simple_toy_MetNet()
    rxns_info_ = Dict()
    for rxn in Chemostat.Utils.rxns(model)
        rxns_info_[rxn] = rand([0,-1, 1])
    end
    test_MetNetUtils_add_met(model, rxns_info_)

    rxns_info_ = Dict()
    for rxn_i in 1:Chemostat.Utils.rxnscount(model)
        rxns_info_[rxn_i] = rand([0,-1, 1])
    end
    test_MetNetUtils_add_met(model, rxns_info_)
end
test_MetNetUtils_add_met()