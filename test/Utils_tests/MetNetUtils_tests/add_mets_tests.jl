function test_add_met(model, rxns_info_)
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

function test_add_met()
    model = Chemostat.Utils.simple_toy_MetNet()
    rxns_info_ = Dict()
    for rxn in Chemostat.Utils.rxns(model)
        rxns_info_[rxn] = rand([0,-1, 1])
    end
    test_add_met(model, rxns_info_)

    rxns_info_ = Dict()
    for rxn_i in 1:Chemostat.Utils.rxnscount(model)
        rxns_info_[rxn_i] = rand([0,-1, 1])
    end
    test_add_met(model, rxns_info_)
end
test_add_met()