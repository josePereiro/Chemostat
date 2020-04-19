"""
    A few useful plot methods
"""
module Plots
    using Chemostat.Utils
    using Distributions
    import Plots: plot, plot!, scatter, scatter!, histogram, histogram!

    include("ep_hr_marginals.jl")
end