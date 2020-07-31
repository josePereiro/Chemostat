module Chemostat

import Requires: @require

include("Utils/Utils.jl")

include("LP/LP.jl")
include("Test/Test.jl")
include("SteadyState/SteadyState.jl")
include("MaxEntEP/MaxEntEP.jl")
include("MaxEntHR/MaxEntHR.jl")
include("Plots/Plots.jl")

end # module