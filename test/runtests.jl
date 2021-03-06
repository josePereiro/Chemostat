using Test
using Random
using Serialization
using Chemostat
using ProgressMeter

@testset "Chemostat.jl" begin

    try
    
        include("Prepare_tests/prepare_tests.jl")
        include("Utils_tests/Utils_tests.jl")
        include("LP_tests/LP_tests.jl")
    catch err
        rethrow(err)
    finally
        # Delete cache
        rm(TEST_CACHE_DIR, force = true, recursive = true)
        @test !isdir(TEST_CACHE_DIR)
    end

end
