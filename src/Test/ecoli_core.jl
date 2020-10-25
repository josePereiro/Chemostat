const ECOLI_MODEL_BIOMASS_IDER = "Biomass_Ecoli_core_w_GAM"

function ecoli_core_model()
    matfile = joinpath(@__DIR__, "..", "..", "data", "test", "ecoli_core_model.mat")
    mat_model = MAT.matread(matfile)["model"]
    return MetNet(mat_model; reshape = true)
end