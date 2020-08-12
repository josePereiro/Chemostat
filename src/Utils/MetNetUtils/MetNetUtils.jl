"""
    A few methods for working with MetNets.
"""

# TODO implement a way to read and write models to csv
# IDEA make a common protocol of linealize all arrays for writing and them 
# reshape for reading
# A = rand(1:10, 5, 6)
# a = vec(A);
# @assert all(reshape(a, 5, 6) .== A)

include("base.jl")
include("defaults.jl")
include("iders.jl")
include("setter.jl")
include("short_loops.jl")
include("getter.jl")
include("queries.jl")
include("read_mat.jl")
include("summary.jl")
include("split_revs.jl")
include("update.jl")
include("add_met.jl")
include("add_rxn.jl")
include("invert_bkwds.jl")
include("make_intake_info_compat.jl")
include("invert_rxn.jl")
include("add_costs.jl")
include("rxn_str.jl")
include("search.jl")
include("interchange_rxn.jl")
include("del_rxn.jl")
include("balance_str.jl")
include("del_blocked.jl")
include("del_met.jl")
include("similar_rxns.jl")
include("compress_metnet.jl")
include("clamp_bounds.jl")