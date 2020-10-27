# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function read_mat(filename::String)
    model_id, model_mat = first(matread(filename))
    return MetNet(model_mat; reshape = true)
end