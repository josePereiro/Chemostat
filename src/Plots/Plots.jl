"""
    A few useful plot methods
"""
module Plots
    import ..Utils: logspace, rxnindex, marginal, rxns, lb, ub, HRout, hists, μ, σ, FBAout, pdf_maxval, av
    import MetabolicEP: EPout, MetNet
    import Distributions: mean, var, pdf, std, truncated, Normal
    import LinearAlgebra: normalize
    import SpecialFunctions: erf
    
    import Plots: plot, plot!, scatter, scatter!, histogram, histogram!

    include("plot_marginal.jl")
end

