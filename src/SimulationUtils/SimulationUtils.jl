
module SimulationUtils

    import Serialization: serialize, deserialize
    import SparseArrays
    import SparseArrays: spzeros, SparseMatrixCSC, SparseVector
    import ..Chemostat.Utils: EPModel, EPout,
                              rxnindex, av, norm1_stoi_err, update_solution!,
                              set_cache_dir, load_cache, save_cache,
                              tagprintln_inmw, temp_cache_file, delete_temp_caches,
                              struct_to_dict, compressed_copy
    
    import ..Chemostat.Test: empty_epout
    import ..Chemostat.MaxEntEP: converge_ep!
    import ..Chemostat.LP: fba
    import StatsBase: mean

    include("cached_simulation.jl")
    include("epoch_converge_ep.jl")
end