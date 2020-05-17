function test_add_get_data()
    model = deserialize(METNET_CACHE_FILE)
    fbaout = deserialize(FBAOUT_CACHE_FILE)
    epout1 = deserialize(EPOUT_ALPHA_FIN_CACHE_FILE)
    epout2 = deserialize(EPOUT_ALPHA_INF_CACHE_FILE)
    
    boundle = Chemostat.Utils.ChstatBoundle()
    ξs = [1, 2]
    βs = [2, 3]
    Chemostat.Utils.add_data!(boundle, ξs[1], :net, model);
    Chemostat.Utils.add_data!(boundle, ξs[1], :fba, fbaout);
    Chemostat.Utils.add_data!(boundle, ξs[1], βs[1], :ep, epout1);
    Chemostat.Utils.add_data!(boundle, ξs[2], βs[2], :ep, epout2);

    @test (Chemostat.Utils.get_data(boundle, ξs[1], :net) === model)
    @test (Chemostat.Utils.get_data(boundle, ξs[1], :fba) === fbaout)
    @test (Chemostat.Utils.get_data(boundle, ξs[1], βs[1], :ep) === epout1)
    @test (Chemostat.Utils.get_data(boundle, ξs[2], βs[2], :ep) === epout2)

end
test_add_get_data()