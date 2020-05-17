function _test_out_getters(model, out)

    av = Chemostat.Utils.av(out)
    va = Chemostat.Utils.va(out)
    μ = Chemostat.Utils.μ(out)
    σ = Chemostat.Utils.σ(out)

    for (i, rxn) in enumerate(Chemostat.Utils.rxns(model))
        @test Chemostat.Utils.av(model, out, rxn) == av[i]
        @test Chemostat.Utils.av(model, out, i) == av[i]
        @test Chemostat.Utils.va(model, out, rxn) == va[i]
        @test Chemostat.Utils.va(model, out, i) == va[i]
        @test Chemostat.Utils.μ(model, out, rxn) == μ[i]
        @test Chemostat.Utils.μ(model, out, i) == μ[i]
        @test Chemostat.Utils.σ(model, out, rxn) == σ[i]
        @test Chemostat.Utils.σ(model, out, i) == σ[i]
    end
    
    rand_idxs = rand(1:length(model.rxns), length(model.rxns))
    @test all(Chemostat.Utils.av(model, out, rand_idxs) .== av[rand_idxs])
    @test all(Chemostat.Utils.va(model, out, rand_idxs) .== va[rand_idxs])
    @test all(Chemostat.Utils.μ(model, out, rand_idxs) .== μ[rand_idxs])
    @test all(Chemostat.Utils.σ(model, out, rand_idxs) .== σ[rand_idxs])

end

#fbaout
function test_fbaout_getters()
    model = deserialize(METNET_CACHE_FILE)
    fbaout = deserialize(FBAOUT_CACHE_FILE)

    @test all(Chemostat.Utils.av(fbaout) .== fbaout.v)
    @test all(Chemostat.Utils.va(fbaout) .== zeros(length(fbaout.v)))
    @test all(Chemostat.Utils.μ(fbaout) .== fbaout.v)
    @test all(Chemostat.Utils.σ(fbaout) .== zeros(length(fbaout.v)))

    _test_out_getters(model, fbaout)

end
test_fbaout_getters()

# epout
function _test_epout_getters(model, epout)

    @test all(Chemostat.Utils.av(epout) .== epout.av)
    @test all(Chemostat.Utils.va(epout) .== epout.va)
    @test all(Chemostat.Utils.μ(epout) .== epout.μ)
    @test all(Chemostat.Utils.σ(epout) .== epout.σ)

    _test_out_getters(model, epout)
end

function test_epout_getters()

    model = deserialize(METNET_CACHE_FILE)
    epout_alpha_fin = deserialize(EPOUT_ALPHA_FIN_CACHE_FILE)
    epout_alpha_inf = deserialize(EPOUT_ALPHA_INF_CACHE_FILE)

    _test_epout_getters(model, epout_alpha_fin)
    _test_epout_getters(model, epout_alpha_inf)
end
test_epout_getters()