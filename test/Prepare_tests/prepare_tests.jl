# const test_cache_dir = "../test_cache"

const TEST_ROOT_DIR = dirname(@__DIR__)
const TEST_CACHE_DIR = joinpath(TEST_ROOT_DIR, "cache")

# cache
if !isdir(TEST_CACHE_DIR)
    mkpath(TEST_CACHE_DIR)
    flush(stdout)
    println("created $TEST_CACHE_DIR")
end

# metent
const METNET_CACHE_FILE = joinpath(TEST_CACHE_DIR, "metnet.jls")
function create_metnet_cache()
    isfile(METNET_CACHE_FILE) && return
    metnet = Chemostat.Utils.toy_model()
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
    epout = Chemostat.MaxEntEP.maxent_ep(metnet, α = Inf)
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
    epout = Chemostat.MaxEntEP.maxent_ep(metnet, α = 1e11)
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
    fbaout = Chemostat.FBA.fba(metnet, "biom")
    serialize(FBAOUT_CACHE_FILE, fbaout)
    flush(stdout)
    println("created $FBAOUT_CACHE_FILE")
    @test true
end
create_fbaout_cache()