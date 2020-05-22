"""
    A few useful plot methods
"""
module Plots

import ..Utils: logspace, rxnindex, marginal, rxns, ChstatBoundle,
    lb, ub, MetNet, EPout, HRout, hists, μ, σ, FBAout, pdf_maxval, av, va, get_data,
    parse_β, parse_ξ, metnet_data_key, AbstractOut,
    av_stoi_err_ep, av_stoi_err_hr, IDER_TYPE, met_rxns
import Distributions: mean, var, pdf, std, truncated, Normal
import LinearAlgebra: normalize
import SpecialFunctions: erf
import Plots: plot, plot!, scatter, scatter!, histogram, histogram!, distinguishable_colors

include("plot_marginal.jl")
# include("plot_xi.jl")
# include("plot_stoi_err_beta.jl")
# include("plot_stoi_err_xi.jl")
# include("plot_mean_stoi_err_beta.jl")
# include("plot_mean_stoi_err_xi.jl")
# include("plot_norm_stoi_err_beta.jl")

end

