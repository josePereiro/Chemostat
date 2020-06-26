"""
    Contains a few Utils for testing
"""
module Test

    import SparseArrays: SparseMatrixCSC
    import ..Utils: MetNet, lb!, add_costs
    import ..LP: fva_preprocess

    include("lineal_model.jl")
    include("toy_model.jl")

end