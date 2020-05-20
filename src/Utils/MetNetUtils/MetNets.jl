fake_metsid(M) = ["M$i" for i in 1:M]
fake_metNames(M) = ["MET $i" for i in 1:M]
fake_metFormulas(M) = ["X" for i in 1:M]


fake_rxnsid(N) = ["r$i" for i in 1:N]
fake_rxnNames(N) = ["RXN $i" for i in 1:N]

# Minimum simple Constructor
function MetNet(S::AbstractMatrix, b::AbstractVector, 
                lb::AbstractVector, ub::AbstractVector, 
                rxns::AbstractVector = fake_rxnsid(size(S,2)), 
                mets::AbstractVector = fake_metsid(size(S,1));
                T = Float64,  
                metNames =  fake_metNames(size(S,1)), 
                metFormulas = fake_metFormulas(size(S,1)), 
                rxnNames = fake_rxnNames(size(S,2)))
    
    M,N = size(S);
    return MetNet{T}(
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
    mat_model = deepcopy(mat_model)

    S = Matrix(mat_model["S"])
    b = collect(vec(mat_model["b"]))
    c = collect(haskey(mat_model, "c") ? vec(mat_model["c"]) : [0.0])
    lb = collect(vec(mat_model["lb"]))
    ub = collect(vec(mat_model["ub"]))
    genes = collect(haskey(mat_model, "genes") ? vec(mat_model["genes"]) : [""])
    rxnGeneMat = collect(haskey(mat_model, "rxnGeneMat") ? Matrix(mat_model["rxnGeneMat"]) : [])
    grRules = collect(haskey(mat_model, "grRules") ? vec(mat_model["grRules"]) : [""])
    mets = collect(haskey(mat_model, "mets") ? vec(mat_model["mets"]) : fake_metsid(size(S, 1)))
    rxns = collect(haskey(mat_model, "rxns") ? vec(mat_model["rxns"]) : fake_rxnsid(size(S, 2)))
    metNames = collect(haskey(mat_model, "metNames") ? vec(mat_model["metNames"]) : [""])
    metFormulas = collect(haskey(mat_model, "metFormulas") ? vec(mat_model["metFormulas"]) : [""])
    rxnNames = collect(haskey(mat_model, "rxnNames") ? vec(mat_model["rxnNames"]) : [""])
    rev = collect(haskey(mat_model, "rev") ? vec(mat_model["rev"]) : [])
    subSystems = collect(haskey(mat_model, "subSystems") ? vec(mat_model["subSystems"]) : [""])
    
    return MetNet{T}(S, b, c, lb, ub, genes, 
        rxnGeneMat, grRules, mets, rxns, metNames, metFormulas, rxnNames, rev, subSystems)
end

"""
    Create a new MetNet from a template but overwriting the fields
    of the template with the passed as kwargs
"""
function MetNet(metnet::MetNet; kwargs...)
    kwargs = Dict(kwargs)
    
    metnet_dict = Dict()
    for field in fieldnames(typeof(metnet))
        metnet_dict[string(field)] = getfield(metnet, field)
    end
    
    for (k, v) in metnet_dict
        sk = Symbol(k)
        if haskey(kwargs, sk)
            metnet_dict[k] = kwargs[sk]
        end
    end
    # TODO include a T_default ?
    T = get(kwargs, :T, Float64)
    return MetNet(metnet_dict, T)
end