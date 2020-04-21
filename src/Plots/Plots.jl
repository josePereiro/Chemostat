"""
    A few useful plot methods
"""
module Plots
    import Chemostat.Utils: rxnindex
    import MetabolicEP: EPout, MetNet
    using Distributions
    import Plots: plot, plot!, scatter, scatter!, histogram, histogram!

    include("ep_hr_marginals.jl")
end