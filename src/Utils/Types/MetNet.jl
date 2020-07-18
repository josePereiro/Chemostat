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
    lb::AbstractVector{T} # fluxes lower bound N elements vector (*)
    ub::AbstractVector{T} # fluxes upper bound N elements vector (*)
    genes::AbstractVector{String} # gene names 
    rxnGeneMat::AbstractMatrix{Float64} # 
    grRules::AbstractVector{String} # gene-reaction rule N elements vector of strings (and / or allowed)
    mets::AbstractVector{String} # metabolites short-name M elements (**)
    rxns::AbstractVector{String} # reactions short-name N elements (*)
    metNames::AbstractVector{String} # metabolites long-names M elements (**)
    metFormulas::AbstractVector{String} # metabolites formula M elements (**)
    rxnNames::AbstractVector{String} # reactions long-names N elements (**)
    rev::AbstractVector{Bool} # reversibility of reactions N elements
    subSystems::AbstractVector # cellular component of fluxes N elements
    intake_info::Dict # Information required for enforcing the Chemostat bounds
end

fake_metsid(M) = ["M$i" for i in 1:M]
fake_metNames(M) = ["MET $i" for i in 1:M]
fake_metFormulas(M) = ["X" for i in 1:M]

fake_rxnsid(N) = ["r$i" for i in 1:N]
fake_rxnNames(N) = ["RXN $i" for i in 1:N]
fake_subSystems(N) = ["Net" for i in 1:N]


# Minimum simple Constructor
function MetNet(S::AbstractMatrix, b::AbstractVector, 
                lb::AbstractVector, ub::AbstractVector, 
                rxns::AbstractVector = fake_rxnsid(size(S,2)), 
                mets::AbstractVector = fake_metsid(size(S,1));
                T = Float64,  
                metNames =  fake_metNames(size(S,1)), 
                metFormulas = fake_metFormulas(size(S,1)), 
                rxnNames = fake_rxnNames(size(S,2)), 
                subSystems = fake_subSystems(size(S,2)),
                intake_info = Dict())
    
    M,N = size(S);
    return MetNet{T}(S, b, zeros(Float64,N), lb, ub, 
        ["NA"], zeros(1,1),
        ["NA"], mets, rxns, 
        metNames, metFormulas, rxnNames, [], 
        ["IN"], intake_info) 
end

function reshape_mat_dict!(mat_dict; S::DataType = String, 
    F::DataType = Float64, N::DataType = Int)
    # String vectors
    for k in ["comps", "metNames", "metFormulas", 
            "rxnFrom", "rxnNames", "genes", "inchis", 
            "grRules", "mets", "rxns"]
        !haskey(mat_dict, k) && continue
        mat_dict[k] = mat_dict[k] |> vec .|> S;
    end

    # Float vectors
    for k in ["b", "lb", "ub", "c"]
        !haskey(mat_dict, k) && continue
        mat_dict[k] = mat_dict[k] |> vec .|> Float64;
    end

    # Integer vectors
    for k in ["metComps"]
        !haskey(mat_dict, k) && continue
        mat_dict[k] = floor.(Int, mat_dict[k] |> vec);
    end

    # S
    if haskey(mat_dict, "S")
        mat_dict["S"] = Matrix{F}(mat_dict["S"])
    end

    # Matrices
    for k in ["rxnGeneMat"]
        !haskey(mat_dict, k) && continue
        mat_dict[k] = mat_dict[k] |> Matrix
    end

    return mat_dict
end

# For COBRA .mat compatible files
function MetNet(mat_model::Dict, T = nothing) 
    mat_model = reshape_mat_dict!(deepcopy(mat_model))

    S = mat_model["S"]
    b = mat_model["b"]
    lb = mat_model["lb"]
    ub = mat_model["ub"]
    mets = get(mat_model, "mets", fake_metsid(size(S, 1)))
    rxns = get(mat_model, "rxns", fake_rxnsid(size(S, 2)))

    return MetNet(S, b, lb, ub, rxns, mets;
            T = isnothing(T) ? eltype(S) : T,  
            metNames =  get(mat_model, "metNames", fake_metNames(size(S,1))), 
            metFormulas = get(mat_model, "metFormulas", fake_metFormulas(size(S,1))), 
            rxnNames = get(mat_model, "rxnNames", fake_rxnNames(size(S,2))), 
            subSystems = get(mat_model, "subSystems", fake_subSystems(size(S,2))),
            intake_info = get(mat_model, "intake_info", Dict()))
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
    T = get(kwargs, :T, nothing)
    return MetNet(metnet_dict, T)
end