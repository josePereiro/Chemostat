
module SimulationUtils

    import Serialization: serialize, deserialize
    import SparseArrays
    import SparseArrays: spzeros, SparseMatrixCSC, SparseVector
    import ..Chemostat.Utils: EPModel, rxnindex, struct_to_dict, compress_dict, 
                              av, norm_abs_stoi_err, update_solution!
    import ..Chemostat.Test: empty_epout
    import ..Chemostat.MaxEntEP: converge_ep!
    import ..Chemostat.LP: fba
    import StatsBase: mean

    include("cache.jl")
    include("print_inmw.jl")
    include("cached_simulation.jl")
    include("epoch_converge_ep.jl")
end