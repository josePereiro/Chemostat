
module SimulationUtils

    import Serialization: serialize, deserialize
    import SparseArrays
    import ProgressMeter: ProgressThresh, Progress, next!, finish!, update!
    import SparseArrays: spzeros, SparseMatrixCSC, SparseVector
    import ..Chemostat.Utils: MetNet, EPModel, EPout,IDER_TYPE,
                              rxnindex, av, norm1_stoi_err, update_solution!,
                              set_cache_dir, load_cache, save_cache,
                              tagprintln_inmw, temp_cache_file, delete_temp_caches,
                              struct_to_dict, compressed_copy, err_str
    
    import ..Chemostat.Test: empty_epout
    import ..Chemostat.MaxEntEP: converge_ep!, maxent_ep
    import ..Chemostat.LP: fba
    import StatsBase: mean

    include("cached_simulation.jl")
    include("epoch_converge_ep.jl")
    include("find_beta.jl")
    
end