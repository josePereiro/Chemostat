module LP

    # TODO MathProgBase is deprecated, update the code!!!
    # Possibly make pull request for MetabolicEP too
    import MathProgBase
    import MathProgBase.HighLevelInterface: linprog
    import Clp: ClpSolver
    import ..Utils: rxnindex, FBAout, MetNet, IDER_TYPE, del_blocked, bounds, fixxing, err_str, IDER_TYPE
    import SparseArrays: issparse
    import ProgressMeter: Progress, update!, finish!, next!
    import Base.Threads: nthreads, threadid, @threads
    import JuMP, GLPK

    include("const.jl")
    include("utils.jl")
    include("fba.jl")
    include("fva.jl")
    include("fva_preprocess.jl")
    include("check_newbounds.jl")
    include("yLP.jl")
    include("projection2D.jl")
    
end

