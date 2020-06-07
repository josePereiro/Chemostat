# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

#=
    The * mark the relevant fields for EP algorithm, 
    and the ** the relevante for other methods in this package,
    the rest are untouched, so they may be contain dummy values sometimes!!!
=#
struct MetNet{T<:Real}
    S::AbstractMatrix{T} # Stoichiometric matrix M x N sparse (*)
    b::AbstractVector{T} # right hand side of equation  S ν = b (*)
    c::AbstractVector{Float64} # reaction index of biomass 
    lb::AbstractVector{T} # fluxes lower bound M elements vector (*)
    ub::AbstractVector{T} # fluxes upper bound M elements vector (*)
    genes::AbstractVector{String} # gene names 
    rxnGeneMat::AbstractMatrix{Float64} # 
    grRules::AbstractVector{String} # gene-reaction rule N elements vector of strings (and / or allowed)
    mets::AbstractVector{String} # metabolites short-name M elements (**)
    rxns::AbstractVector{String} # reactions short-name N elements (*)
    metNames::AbstractVector{String} # metabolites long-names M elements (**)
    metFormulas::AbstractVector{String} # metabolites formula M elements (**)
    rxnNames::AbstractVector{String} # reactions long-names N elements (**)
    rev::AbstractVector{Bool} # reversibility of reactions N elements
    subSystems::AbstractVector{String} # cellular component of fluxes N elements
    intake_info::Dict # Information required for enforcing the Chemostat bounds
end

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
    return MetNet{T}(S, b, zeros(Float64,N), lb, ub, 
        ["NA"], SparseMatrixCSC{Float64,Int64}(zeros(1,1)),
        ["NA"], mets, rxns, metNames, metFormulas, rxnNames, (lb .< 0.0) .& (ub .> 0.0), 
        ["IN"], Dict()) 
end

# For COBRA .mat compatible files
function MetNet(mat_model::Dict, T = Float64) 
    mat_model = deepcopy(mat_model)

    S = Matrix(mat_model["S"])
    b = collect(vec(mat_model["b"]))
    c = get(mat_model, "c", [0.0]) |> vec |> collect
    lb = mat_model["lb"] |> vec |> collect
    ub = mat_model["ub"] |> vec |> collect
    genes = get(mat_model, "genes", [""]) |> vec |> collect
    rxnGeneMat = get(mat_model, "rxnGeneMat", []) |> Matrix |> collect
    grRules = get(mat_model, "grRules", [""]) |> vec |> collect
    mets = get(mat_model, "mets", fake_metsid(size(S, 1))) |> vec |> collect
    rxns = get(mat_model, "rxns", fake_rxnsid(size(S, 2))) |> vec |> collect
    metNames = get(mat_model, "metNames", [""]) |> vec |> collect
    metFormulas = get(mat_model, "metFormulas", [""]) |> vec |> collect
    rxnNames = get(mat_model, "rxnNames", [""]) |> vec |> collect
    rev = get(mat_model, "rev", []) |> vec |> collect
    subSystems = get(mat_model, "subSystems", [""]) |> vec |> collect
    intake_info = get(mat_model, "intake_info", Dict())

    return MetNet{T}(S, b, c, lb, ub, genes, 
        rxnGeneMat, grRules, mets, rxns, metNames, 
        metFormulas, rxnNames, rev, subSystems, intake_info)
end

"""
    Create a new MetNet from a template but overwriting the fields
    of the template with the given as kwargs
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