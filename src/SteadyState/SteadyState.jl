module SteadyState

import Chemostat.Utils: rxnindex, lb!, ub!, rxns
import MetabolicEP: MetNet

include("stst_bound.jl")
include("apply_bound.jl")

end