# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

"""
    Common utilities
"""
module Utils
    
    import MathProgBase.HighLevelInterface: linprog
    import Clp: ClpSolver
    import SparseArrays: SparseMatrixCSC, spzeros, findnz, sparsevec, sparse
    import StatsBase: AbstractHistogram, Histogram, fit
    import Distributions: Normal, Truncated, mean, var, pdf
    import LinearAlgebra: normalize, qr, diag
    import Base: isequal, size, ==
    import JSON: json, parse
    
    include("General/General.jl")
    include("Types/Types.jl")
    # include("MetabolicEPUtils/MetabolicEPUtils.jl") # Deprecated
    include("MetNetUtils/MetNetUtils.jl")
    include("OutUtils/OutUtils.jl")
    include("TestUtils/TetsUtils.jl")
    include("BoundleUtils/BoundleUtils.jl")
    include("MaxEntEPUtils/MaxEntUtils.jl")
end
