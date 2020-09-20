function test_get_set_index()
    model = deserialize(METNET_CACHE_FILE)
    fbaout = deserialize(FBAOUT_CACHE_FILE)
    epout1 = deserialize(EPOUT_ALPHA_FIN_CACHE_FILE)
    epout2 = deserialize(EPOUT_ALPHA_INF_CACHE_FILE)
    
    bundle = Chemostat.Utils.ChstatBundle()
    ξs = [1, 2]
    βs = [2, 3]
    bundle[ξs[1], :net] = model;
    bundle[ξs[1], :fba] = fbaout;
    bundle[ξs[1], βs[1], :ep] = epout1;
    bundle[ξs[2], βs[2], :ep] = epout2;

    @test (bundle[ξs[1], :net] === model)
    @test (bundle[ξs[1], :fba] === fbaout)
    @test (bundle[ξs[1], βs[1], :ep] === epout1)
    @test (bundle[ξs[2], βs[2], :ep] === epout2)
end
test_get_set_index()