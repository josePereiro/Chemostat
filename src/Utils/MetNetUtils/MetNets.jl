fake_metsid(M) = ["M$i" for i in 1:M]
fake_metNames(M) = ["MET $i" for i in 1:M]
fake_metFormulas(M) = ["X" for i in 1:M]


fake_rxnsid(N) = ["r$i" for i in 1:N]
fake_rxnNames(N) = ["RXN $i" for i in 1:N]

# Minimum Constructor
function MetNet(S, b, lb, ub, 
    xns = fake_rxnsid(size(S,2)), 
    mets = fake_metsid(size(S,1));
    T = Float64,  
    metNames =  fake_metNames(size(S,1)), 
    metFormulas = fake_metFormulas(size(S,1)), 
    rxnNames = fake_rxnNames(size(S,2)))
    
    M,N = size(S);
    return MetNet{T}(
        N, # number of fluxes
        M, # number of metabolites
        S, # Stoichiometric matrix M x N sparse
        b, # right hand side of equation  S ν = b 
        zeros(Float64,N), # reaction index of biomass
        lb, # fluxes lower bound N elements vector
        ub, # fluxes upper bound N elements vector
        ["NA"], # gene names
        SparseMatrixCSC{Float64,Int64}(zeros(1,1)), # rxnGeneMat
        ["NA"], # gene-reaction rule N elements vector of strings (and / or allowed)
        mets, # metabolites short-name M elements 
        rxns, # reactions short-name N elements
        metNames, # metabolites long-names M elements
        metFormulas, # metabolites formula M elements
        rxnNames, # reactions long-name N elements
        (lb .< 0.0) .& (ub .> 0.0), # reversibility of reactions N elements
        ["IN"]) # cellular component of fluxes N elements
end

# For COBRA .mat compatible files
function MetNet(mat_model::Dict, T = Float64) 

    S = Matrix(mat_model["S"])
    M, N = size(S)
    b = reshape(mat_model["b"], :, 1)
    c = haskey(mat_model, "c") ? reshape(mat_model["c"], :, 1) : [0.0]
    lb = reshape(mat_model["lb"], :, 1)
    ub = reshape(mat_model["ub"], :, 1)
    genes = haskey(mat_model, "genes") ? String.(reshape(mat_model["genes"], :, 1)) : [""]
    rxnGeneMat = haskey(mat_model, "rxnGeneMat") ? Matrix(mat_model["rxnGeneMat"]) : []
    grRules = haskey(mat_model, "grRules") ? String.(reshape(mat_model["grRules"], :, 1)) : [""]
    mets = haskey(mat_model, "mets") ? String.(reshape(mat_model["mets"], :, 1)) : ["M$i" for i in 1:M]
    rxns = haskey(mat_model, "rxns") ? String.(reshape(mat_model["rxns"], :, 1)) : ["r$i" for i in 1:N]
    metNames = haskey(mat_model, "metNames") ? String.(reshape(mat_model["metNames"], :, 1)) : [""]
    metFormulas = haskey(mat_model, "metFormulas") ? String.(reshape(mat_model["metFormulas"], :, 1)) : [""]
    rxnNames = haskey(mat_model, "rxnNames") ? String.(reshape(mat_model["rxnNames"], :, 1)) : [""]
    rev = haskey(mat_model, "rev") ? Bool.(reshape(mat_model["rev"], :, 1)) : (lb .< 0.0) .& (ub .> 0.0)
    subSystems = haskey(mat_model, "subSystems") ? String.(reshape(mat_model["subSystems"], :, 1)) : [""]
    
    return MetNet{T}(N, M, S, b, c, lb, ub, genes, 
        rxnGeneMat, grRules, mets, rxns, metNames, metFormulas, rxnNames, rev, subSystems)
end