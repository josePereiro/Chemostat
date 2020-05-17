using Chemostat
using Test
using Random
using Serialization

@testset "Chemostat.jl" begin
    
    include("Prepare_tests/prepare_tests.jl")
    include("Utils_tests/Utils_tests.jl")

end
