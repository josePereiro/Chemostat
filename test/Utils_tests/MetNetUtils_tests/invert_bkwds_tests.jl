function tests_invert_bkwds()
    model = Chemostat.Utils.simple_toy_MetNet()
    bkwd_idxs = Chemostat.Utils.bkwds(model)
    new_model = Chemostat.Utils.invert_bkwds(model)

    for rxn_i in bkwd_idxs
        @test all(new_model.S[:,rxn_i] .== -model.S[:,rxn_i])
        @test new_model.lb[rxn_i] == abs(model.ub[rxn_i])
        @test new_model.ub[rxn_i] == abs(model.lb[rxn_i])
        @test new_model.rxns[rxn_i] == model.rxns[rxn_i] * Chemostat.Utils.bkwd_prefix
        @test new_model.rxnNames[rxn_i] == model.rxnNames[rxn_i]
    end
end
tests_invert_bkwds()