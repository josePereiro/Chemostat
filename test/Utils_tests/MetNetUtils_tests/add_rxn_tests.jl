function test_add_rxn(model, mets_info_)
    rxnid_ = "new_rxn"
    rxnName_ = "Metallica!!!"
    c_ = 123
    lb_ = -12
    ub_ = 345
    new_model = Chemostat.Utils.add_rxn(model, rxnid_, 
            mets = mets_info_, rxnName = rxnName_,
            c = c_, lb = lb_, ub = ub_)

    @test new_model.rxns[end] == rxnid_
    @test new_model.rxnNames[end] == rxnName_
    @test new_model.c[end] == c_
    @test new_model.lb[end] == lb_
    @test new_model.ub[end] == ub_
    @test_throws ErrorException Chemostat.Utils.add_rxn(model, model.rxns[1])

    for (met, stoi_coe) in mets_info_
        met_indx = Chemostat.Utils.metindex(new_model, met)
        @test new_model.S[met_indx, end] == stoi_coe
    end
end
function test_add_rxn()
    model = Chemostat.Utils.simple_toy_MetNet()

    mets_info_ = Dict()
    for met in Chemostat.Utils.mets(model)
        mets_info_[met] = rand([0,-1, 1])
    end
    test_add_rxn(model, mets_info_)

    mets_info_ = Dict()
    for met_i in 1:Chemostat.Utils.metscount(model)
        mets_info_[met_i] = rand([0,-1, 1])
    end
    test_add_rxn(model, mets_info_)
end
test_add_rxn()