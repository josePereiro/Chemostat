"""
    Contains a few Utils for testing
"""
module Test

    import SparseArrays: SparseMatrixCSC
    import ..Utils: MetNet, lb!, EPout, EPFields, expanded_model, 
                    Met, Rxn, set_met!, set_rxn!
    import ..LP: fva_preprocess
    import MAT

    include("lineal_model.jl")
    include("toy_model.jl")
    include("ecoli_core.jl")
    include("test_types.jl")

end