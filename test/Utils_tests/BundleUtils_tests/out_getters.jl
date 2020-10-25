function test_out_getters()
    bundle = deserialize(BUNDLE_CACHE_FILE)

    getters = [Chemostat.Utils.av, Chemostat.Utils.va, 
        Chemostat.Utils.μ, Chemostat.Utils.σ]
    
    ξs = bundle.ξs
    βs = bundle.βs 

    for getter in getters

        for (ξi, ξ) in enumerate(ξs)
            model = bundle[ξ, :net]
            fbaout = bundle[ ξ, :fba]
            
            # Indexing using ξ, data_key
            @test all(getter(bundle, ξ, :fba) .== getter(fbaout))
            
            # Indexing using ξ, data_key, iders
            rand_idxs = rand(1:size(model, 2), size(model, 2))
            @test all(getter(bundle, ξ, :fba, rand_idxs) .== 
                            getter(model, fbaout, rand_idxs))
            @test all(getter(bundle, ξ, :fba, model.rxns[rand_idxs]) .== 
                            getter(model, fbaout, model.rxns[rand_idxs]))

            # Indexing using ξs, data_key, ider
            @test all(getter(bundle, ξs, :fba, rand_idxs[1])[ξi] .== 
                            getter(model, fbaout, rand_idxs[1]))

            for (βi, β) in enumerate(βs)
                epout = bundle[ξ, β, :ep]

                # Indexing using ξ, β, data_key
                @test all(getter(bundle, ξ, β, :ep) .== getter(epout))

                # Indexing using ξ, β, data_key, iders
                @test all(getter(bundle, ξ, β, :ep, rand_idxs) .== 
                            getter(model, epout, rand_idxs))
                @test all(getter(bundle, ξ, β, :ep, model.rxns[rand_idxs]) .== 
                            getter(model, epout, model.rxns[rand_idxs]))

                # Indexing using ξs, β, data_key, ider
                @test all(getter(bundle, ξs, β, :ep, rand_idxs[1])[ξi] .== 
                            getter(model, epout, rand_idxs[1]))
                @test all(getter(bundle, ξ, βs, :ep, rand_idxs[1])[βi] .== 
                            getter(model, epout, rand_idxs[1]))

            end
        end
    end
end
test_out_getters()