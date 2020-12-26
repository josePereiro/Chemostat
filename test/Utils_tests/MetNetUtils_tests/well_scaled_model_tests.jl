function test_well_scaled_model()
    scale_factors = 10.0.^(0.005:0.001:0.1)
    orig_model = Chemostat.Test.ecoli_core_model()
    obj_ider = Chemostat.Test.ECOLI_MODEL_BIOMASS_IDER

    prog = Progress(length(scale_factors))
    for scale_factor in scale_factors

        scl_model = Chemostat.Utils.well_scaled_model(orig_model, scale_factor; 
            verbose = true)

        orig_fbaout = Chemostat.LP.fba(orig_model, obj_ider)
        scl_fbaout = Chemostat.LP.fba(scl_model, obj_ider)

        for orig_rxn in orig_model.rxns
            scl_rxn = Chemostat.Utils.rxnindex(scl_model, orig_rxn)
            orig_av = Chemostat.Utils.av(orig_model, orig_fbaout, orig_rxn)
            scl_av = Chemostat.Utils.av(scl_model, scl_fbaout, scl_rxn)

            @test isapprox(orig_av, scl_av; atol = 1e-8)
        end

        next!(prog; showvalues = [
                ("scale_factor        ", scale_factor),
                (": ----------------- ", "ORIGINAL MODEL"),
                ("model size:         ", size(orig_model)),
                ("nzabs_range:        ", Chemostat.Utils.nzabs_range(orig_model.S)),
                ("obj_val:            ", Chemostat.Utils.av(orig_model, orig_fbaout, obj_ider)),
                (": ----------------- ", "SCALED MODEL"),
                ("model size:         ", size(scl_model)),
                ("nzabs_range:        ", Chemostat.Utils.nzabs_range(scl_model.S)),
                ("obj_val:            ", Chemostat.Utils.av(scl_model, scl_fbaout, obj_ider))
            ]
        )
    end
    finish!(prog)
end
test_well_scaled_model()