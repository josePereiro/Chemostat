"""
    A few useful plot methods
"""
module Plots

import ..Utils: logspace, rxnindex, marginal, rxns, ChstatBoundle,
    lb, ub, HRout, hists, μ, σ, FBAout, pdf_maxval, av, va, get_data,
    get_metnet, get_epout, get_fbaout, get_hrout, parse_β, parse_ξ,
    av_stoi_err_ep, av_stoi_err_hr, IDER_TYPE
import MetabolicEP: EPout, MetNet
import Distributions: mean, var, pdf, std, truncated, Normal
import LinearAlgebra: normalize
import SpecialFunctions: erf
import Plots: plot, plot!, scatter, scatter!, histogram, histogram!, distinguishable_colors

include("plot_marginal.jl")
include("plot_xi.jl")
include("plot_stoi_err_beta.jl")
include("plot_stoi_err_xi.jl")
include("plot_mean_stoi_err_beta.jl")
include("plot_mean_stoi_err_xi.jl")

end

