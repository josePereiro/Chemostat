#TODO reference Cossios MaxEnt paper
module MaxEnt

import ..Utils: rxnindex, hrsamples, HRout, eachsample, has_samples, IDER_TYPE
import MetabolicEP: EPout, MetNet, EPMat, EPMatT0
import MetabolicEP.HitAndRun: hrsample
import SparseArrays: spzeros
import LinearAlgebra: diag, Diagonal, mul!
import Distributions: Truncated, Normal, mean, var

include("maxent_hrsample.jl")
include("maxent_metabolicEP.jl")
include("maxent_metabolicEPT0.jl")

end
