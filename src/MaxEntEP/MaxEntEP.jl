
#= The EP implementation is derived from:
Alfredo Braunstein, Anna Paola Muntoni, and Andrea Pagnani, “An Analytic Approximation of the Feasible 
Space of Metabolic Networks,” Nature Communications 8 (April 6, 2017): https://doi.org/10.1038/ncomms14915.
Package metabolicEP (https://github.com/anna-pa-m/Metabolic-EP) 
It was modified just for implement MaxEnt as described in the Appendix of:
Jorge Fernandez-de-Cossio-Diaz and Roberto Mulet, “Maximum Entropy and Population Heterogeneity in 
Continuous Cell Cultures.” PLOS Computational Biology 15, no. 2 (February 27, 2019): e1006823. 
https://doi.org/10.1371/journal.pcbi.1006823.
=#

module MaxEntEP

import SparseArrays: SparseMatrixCSC
import ExtractMacro: @extract
import LinearAlgebra: diag, rmul!, inv!, Hermitian, cholesky!, Diagonal, mul!, diagind
import Printf: @printf
import SpecialFunctions: erf
import Distributions: Truncated, Normal, mean, var
import ProgressMeter: Progress, update!, finish!, ProgressThresh

import ..Utils: AbstractEPMat, EPAlg, EPFields, EPMat, EPMatT0, EPout, MetNet

include("compute_mom5d.jl")
include("epconverge.jl")
include("eponesweep.jl")
include("eponesweepT0.jl")
include("fast_maxent_ep.jl")
include("matchmom.jl")
include("maxent_ep.jl")
include("newav.jl")
include("newμs.jl")
include("prepare_beta_v.jl")
include("prepareinput.jl")
include("scaleepfield.jl")
include("utils.jl")
# include("Q_sigma.jl")

end