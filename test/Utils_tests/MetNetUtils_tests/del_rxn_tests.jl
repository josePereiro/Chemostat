function test_del_rxn()
    model = deserialize(METNET_CACHE_FILE)

    M, N = size(model)
    model.S .= rand(M, N)
    model.c .= rand(N)
    model.lb .= -rand(N)
    model.ub .= rand(N)

    for rxn in model.rxns
        Chemostat.Utils.del_rxn!(model, rxn);
    end

    @test iszero(model.S) 
    @test iszero(model.lb) 
    @test iszero(model.ub) 
    @test iszero(model.c) 
    @test all(model.rxns .== Chemostat.Utils.EMPTY_SPOT) 
end
test_del_rxn()