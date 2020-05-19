# TODO test the rests of the fields
function _test_del_rxn(model, ider)
    test_idx2 = Chemostat.Utils.rxnindex(model, ider)
    test_idx1 = test_idx2 + 1
    if test_idx2 == size(model, 2)
        test_idx2 = test_idx2 - 1
        test_idx1 = test_idx2
    end
    
    del_model = Chemostat.Utils.del_rxn(model, ider);
    @test all(model.S[:, test_idx1] .== del_model.S[:, test_idx2])
    @test model.lb[test_idx1] == del_model.lb[test_idx2]
    @test model.ub[test_idx1] == del_model.ub[test_idx2]
    @test model.rxns[test_idx1] == del_model.rxns[test_idx2]
    
end
function test_del_rxn()
    model = deserialize(METNET_CACHE_FILE)

    for (rxni, rxn) in enumerate(model.rxns) 
        _test_del_rxn(model, rxn)
        _test_del_rxn(model, rxni)
    end
end
test_del_rxn()