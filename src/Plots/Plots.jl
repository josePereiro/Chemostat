"""
    A few useful plot methods
"""
module Plots

import ..Utils: logspace, rxnindex, marginal, rxns, ChstatBoundle,
    lb, ub, MetNet, EPout, HRout, hists, μ, σ, FBAout, pdf_maxval, av, va,
    parse_β, parse_ξ, metnet_data_key, AbstractOut, IDER_TYPE, met_rxns, collect_data
import Distributions: mean, var, pdf, std, truncated, Normal
import LinearAlgebra: normalize
import Plots: plot, plot!, scatter, scatter!, histogram, histogram!, distinguishable_colors

include("plot_marginal.jl")
# include("plot_xi.jl")
# include("plot_stoi_err_beta.jl")
# include("plot_stoi_err_xi.jl")
# include("plot_mean_stoi_err_beta.jl")
# include("plot_mean_stoi_err_xi.jl")
# include("plot_norm_stoi_err_beta.jl")

end

