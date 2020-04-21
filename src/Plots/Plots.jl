"""
    A few useful plot methods
"""
module Plots
    import Chemostat.Utils: rxnindex, trunc_normal
    import MetabolicEP: EPout, MetNet
    import Distributions: mean, var, pdf
    import Plots: plot, plot!, scatter, scatter!, histogram, histogram!

    include("ep_hr_marginals.jl")
end