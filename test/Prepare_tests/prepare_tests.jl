
const TEST_ROOT_DIR = dirname(@__DIR__)
const TEST_CACHE_DIR = joinpath(TEST_ROOT_DIR, "cache")
rm(TEST_CACHE_DIR, force = true, recursive = true)
@test !isdir(TEST_CACHE_DIR)

# cache
if !isdir(TEST_CACHE_DIR)
    mkpath(TEST_CACHE_DIR)
    @test isdir(TEST_CACHE_DIR)
    flush(stdout)
    println("created $TEST_CACHE_DIR")
end

# metent
const METNET_CACHE_FILE = joinpath(TEST_CACHE_DIR, "metnet.jls")
function create_metnet_cache()
    isfile(METNET_CACHE_FILE) && return
    metnet = Chemostat.Test.toy_model()
    serialize(METNET_CACHE_FILE, metnet)
    flush(stdout)
    println("created $METNET_CACHE_FILE")
    @test true
end
create_metnet_cache()

#epout
const EPOUT_ALPHA_INF_CACHE_FILE = joinpath(TEST_CACHE_DIR, "epout_alpha_inf.jls")
function create_epout_alpha_inf_cache()
    isfile(EPOUT_ALPHA_INF_CACHE_FILE) && return
    metnet = deserialize(METNET_CACHE_FILE)
    epout = Chemostat.MaxEntEP.maxent_ep(metnet, alpha = Inf)
    serialize(EPOUT_ALPHA_INF_CACHE_FILE, epout)
    flush(stdout)
    println("created $EPOUT_ALPHA_INF_CACHE_FILE")
    @test true
end
create_epout_alpha_inf_cache()

const EPOUT_ALPHA_FIN_CACHE_FILE = joinpath(TEST_CACHE_DIR, "epout_alpha_fin.jls")
function create_epout_alpha_fin_cache()
    isfile(EPOUT_ALPHA_FIN_CACHE_FILE) && return
    metnet = deserialize(METNET_CACHE_FILE)
    epout = Chemostat.MaxEntEP.maxent_ep(metnet, alpha = 1e11)
    serialize(EPOUT_ALPHA_FIN_CACHE_FILE, epout)
    flush(stdout)
    println("created $EPOUT_ALPHA_FIN_CACHE_FILE")
    @test true
end
create_epout_alpha_fin_cache()

#fbaout
const FBAOUT_CACHE_FILE = joinpath(TEST_CACHE_DIR, "fbaout.jls")
function create_fbaout_cache()
    isfile(FBAOUT_CACHE_FILE) && return
    metnet = deserialize(METNET_CACHE_FILE)
    fbaout = Chemostat.LP.fba(metnet, "biom")
    serialize(FBAOUT_CACHE_FILE, fbaout)
    flush(stdout)
    println("created $FBAOUT_CACHE_FILE")
    @test true
end
create_fbaout_cache()

#Bundle
const BUNDLE_CACHE_FILE = joinpath(TEST_CACHE_DIR, "chstat_bundle.jls")
function create_bundle_cache()
    # This is a whole workflow

    βs = [0.0; Chemostat.Utils.logspace(-2,5, 10)];
    ξs = Chemostat.Utils.logspace(-1,1, 10);

    # this only works for the toymodel
    intake_info = Dict(
        # Open intakes
        "gt" => Dict("ub" => 100.0, "c" => 10.0),
    );
    obj_ider = "biom"; 

    bundle = Chemostat.Utils.ChstatBundle();

    verbose_ = false
    for (ξi, ξ) in enumerate(ξs)
        
        # Model
        model = deserialize(METNET_CACHE_FILE)
        Chemostat.SteadyState.apply_bound!(model, ξ, intake_info)
        model = Chemostat.LP.fva_preprocess(model, verbose = verbose_)
        
        # FBA
        fbaout = Chemostat.LP.fba(model, obj_ider)
        
        # Add to bundle
        bundle[ξ, :net] =  model
        bundle[ξ, :fba] = fbaout
        
        # MaxEnt-EP
        βv = zeros(size(model, 2))
        obj_idx = Chemostat.Utils.rxnindex(model, obj_ider)

        for (βi, β) in enumerate(βs)

            print("xi: [$(ξi)/ $(length(ξs))] beta: [$(βi)/ $(length(βs))] \r");flush(stdout)
                
            βv[obj_idx] = β
            epout = Chemostat.MaxEntEP.maxent_ep(model; alpha = 1e11, beta_vec = βv, verbose = verbose_)
            bundle[ξ, β, :ep] = epout
            @test true
        end
    end
    println("Done                              ");

    serialize(BUNDLE_CACHE_FILE, bundle)
    flush(stdout)
    println("created $BUNDLE_CACHE_FILE")
    @test true

end
create_bundle_cache()