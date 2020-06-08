function test_get_set_index()
    model = deserialize(METNET_CACHE_FILE)
    fbaout = deserialize(FBAOUT_CACHE_FILE)
    epout1 = deserialize(EPOUT_ALPHA_FIN_CACHE_FILE)
    epout2 = deserialize(EPOUT_ALPHA_INF_CACHE_FILE)
    
    boundle = Chemostat.Utils.ChstatBoundle()
    ξs = [1, 2]
    βs = [2, 3]
    boundle[ξs[1], :net] = model;
    boundle[ξs[1], :fba] = fbaout;
    boundle[ξs[1], βs[1], :ep] = epout1;
    boundle[ξs[2], βs[2], :ep] = epout2;

    @test (boundle[ξs[1], :net] === model)
    @test (boundle[ξs[1], :fba] === fbaout)
    @test (boundle[ξs[1], βs[1], :ep] === epout1)
    @test (boundle[ξs[2], βs[2], :ep] === epout2)
end
test_get_set_index()