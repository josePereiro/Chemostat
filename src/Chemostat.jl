module Chemostat

import Requires: @require

include("Utils/Utils.jl")

include("LP/LP.jl")
include("Test/Test.jl")
include("SteadyState/SteadyState.jl")
include("MaxEntEP/MaxEntEP.jl")
include("MaxEntHR/MaxEntHR.jl")
include("SimulationUtils/SimulationUtils.jl")

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        include("Plots/Plots.jl")
    end
end

end # module