"""
    A few useful plot methods
"""
module Plots
    import ..Utils: logspace, rxnindex, marginal, rxns, lb, ub, HRout, hists, μ, σ, FBAout, pdf_maxval
    import MetabolicEP: EPout, MetNet
    import Distributions: mean, var, pdf
    import LinearAlgebra: normalize
    import Plots: plot, plot!, scatter, scatter!, histogram, histogram!

    include("plot_marginal.jl")
    include("normal_xs_for_plotting.jl")
end

