function test_fva_consistency()

    # model
    model = Chemostat.Test.ecoli_core_model();
    m, n = size(model)

    # parameters
    niter = 100
    robj = Chemostat.Test.ECOLI_MODEL_BIOMASS_IDER
    
    # Testing fva consistency checking and not checking an obj value
    for (testname, check_obj) in [("No-Checking Obj, FVA  ", nothing), 
                                  ("Checking Obj, FVA  ", robj)]

        prog = Progress(niter; desc = testname)
        ref_lbs, ref_ubs = Chemostat.LP.fva(model; check_obj, verbose = false)
        for i in 1:niter
            lbs, ubs = Chemostat.LP.fva(model; check_obj, verbose = false)
            @test all(isapprox.(lbs, ref_lbs))
            @test all(isapprox.(ubs, ref_ubs))
            next!(prog)
        end
        finish!(prog)
    end

end
test_fva_consistency()