@testset "MetNetUtils.jl" begin

    include("iders_tests.jl")
    include("getters_tests.jl")
    include("set_mets_tests.jl")
    include("set_rxn_tests.jl")
    include("invert_bkwds_tests.jl")
    include("del_rxn_tests.jl")
    include("well_scaled_model_tests.jl")

end