#= The HitAndRun implementation is derived from:
Package metabolicEP (https://github.com/anna-pa-m/Metabolic-EP) 
It was modified just for implement MaxEnt as described in:
Jorge Fernandez-de-Cossio-Diaz and Roberto Mulet, “Maximum Entropy and Population Heterogeneity in 
Continuous Cell Cultures.” PLOS Computational Biology 15, no. 2 (February 27, 2019): e1006823. 
https://doi.org/10.1371/journal.pcbi.1006823.
=#

module MaxEntHR

import LinearAlgebra: norm
import ..Utils: rxnindex, HRAlg, HRout, eachsample, has_samples, IDER_TYPE, MetNet

include("hrsample.jl")
include("maxent_hr.jl")

end #end module HitAndRun
