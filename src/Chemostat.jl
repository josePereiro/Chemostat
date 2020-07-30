module Chemostat

import Requires: @require

include("Utils/Utils.jl")

include("LP/LP.jl")
include("Test/Test.jl")
include("SteadyState/SteadyState.jl")
include("MaxEntEP/MaxEntEP.jl")
include("MaxEntHR/MaxEntHR.jl")
@require Plots="4076af6c-e467-56ae-b986-b466b2749572" begin
    include("Plots/Plots.jl")
end

end # module