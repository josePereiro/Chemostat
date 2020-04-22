#TODO reference Cossios MaxEnt paper
module MaxEnt

import ..Utils: rxnindex
import MetabolicEP: EPout, MetNet, EPMatT0, AbstractEPMat
import MetabolicEP.HitAndRun: hrsample
import SparseArrays: spzeros
import LinearAlgebra: diag
import Distributions: Truncated, Normal, mean, var

include("maxent_hrsample.jl")
include("maxent_metabolicEP.jl")

end
