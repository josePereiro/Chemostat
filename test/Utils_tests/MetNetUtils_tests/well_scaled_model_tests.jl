function test1_well_scaled_model()
    scale_factors = 10.0.^(0.005:0.001:0.1)
    orig_model = Chemostat.Test.ecoli_core_model()
    orig_min_order, orig_max_order = Chemostat.Utils.nzabs_range(orig_model.S)
    obj_ider = Chemostat.Test.ECOLI_MODEL_BIOMASS_IDER

    prog = Progress(length(scale_factors); desc = "Testing well scaling ...")
    for scale_factor in scale_factors

        scl_model = Chemostat.Utils.well_scaled_model(orig_model, scale_factor; 
            verbose = true)
        scl_min_order, scl_max_order = Chemostat.Utils.nzabs_range(scl_model.S)

        @test scl_min_order >= orig_min_order
        @test scl_max_order <= orig_max_order

        orig_fbaout = Chemostat.LP.fba(orig_model, obj_ider)
        scl_fbaout = Chemostat.LP.fba(scl_model, obj_ider)

        for rxn in orig_model.rxns
            orig_av = Chemostat.Utils.av(orig_model, orig_fbaout, rxn)
            scl_av = Chemostat.Utils.av(scl_model, scl_fbaout, rxn)

            @test isapprox(orig_av, scl_av; atol = 1e-8)
        end

        next!(prog; showvalues = [
                ("scale_factor        ", scale_factor),
                ("                    ", ""),
                (" ORIGINAL MODEL     ", ""),
                ("model size:         ", size(orig_model)),
                ("nzabs_range:        ", (orig_min_order, orig_max_order)),
                ("obj_val:            ", Chemostat.Utils.av(orig_model, orig_fbaout, obj_ider)),
                ("                    ", ""),
                (" SCALED MODEL       ", ""),
                ("model size:         ", size(scl_model)),
                ("nzabs_range:        ", (scl_min_order, scl_max_order)),
                ("obj_val:            ", Chemostat.Utils.av(scl_model, scl_fbaout, obj_ider))
            ]
        )
    end
    finish!(prog)
end
test1_well_scaled_model()

## ------------------------------------------------------------------
# This test just test the same obj_fun value but test all reactions as obj_fun 
function test2_well_scaled_model()

    scale_factors = 10.0.^(0.005:0.005:0.1)
    orig_model = Chemostat.Test.ecoli_core_model()
    # Open all exchanges
    exchs = Chemostat.Utils.exchanges(orig_model)
    Chemostat.Utils.lb!.([orig_model], exchs, -1000.0)

    orig_min_order, orig_max_order = Chemostat.Utils.nzabs_range(orig_model.S)

    N = length(scale_factors) * size(orig_model, 2)
    prog = Progress(N; desc = "Testing well scaling ...")
    for scale_factor in scale_factors

        scl_model = Chemostat.Utils.well_scaled_model(orig_model, scale_factor; 
            verbose = true)
        scl_min_order, scl_max_order = Chemostat.Utils.nzabs_range(scl_model.S)

        @test scl_min_order >= orig_min_order
        @test scl_max_order <= orig_max_order

        for obj_ider in orig_model.rxns
        
            orig_fbaout = Chemostat.LP.fba(orig_model, obj_ider)
            scl_fbaout = Chemostat.LP.fba(scl_model, obj_ider)

            orig_av = Chemostat.Utils.av(orig_model, orig_fbaout, obj_ider)
            scl_av = Chemostat.Utils.av(scl_model, scl_fbaout, obj_ider)

            @test isapprox(orig_av, scl_av; atol = 1e-8)

            next!(prog; showvalues = [
                    ("scale_factor        ", scale_factor),
                    ("obj_ider            ", obj_ider),
                    ("                    ", ""),
                    (" ORIGINAL MODEL     ", ""),
                    ("model size:         ", size(orig_model)),
                    ("nzabs_range:        ", (orig_min_order, orig_max_order)),
                    ("obj_val:            ", Chemostat.Utils.av(orig_model, orig_fbaout, obj_ider)),
                    ("                    ", ""),
                    (" SCALED MODEL       ", ""),
                    ("model size:         ", size(scl_model)),
                    ("nzabs_range:        ", (scl_min_order, scl_max_order)),
                    ("obj_val:            ", Chemostat.Utils.av(scl_model, scl_fbaout, obj_ider))
                ]
            )
        end
    end
    finish!(prog)
end
test2_well_scaled_model()