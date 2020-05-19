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
end