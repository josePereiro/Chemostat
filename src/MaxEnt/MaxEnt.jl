#TODO reference Cossios MaxEnt paper
module MaxEnt

import Chemostat.Utils: rxnindex
import MetabolicEP: EPout, MetNet
import MetabolicEP.HitAndRun: hrsample

using Distributions

include("maxent_hrsample.jl")

end
