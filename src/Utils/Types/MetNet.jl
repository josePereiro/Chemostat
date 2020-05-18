# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

#=
    The * mark the relevant fields for EP algorithm, 
    and the ** the relevante for other methods in this package,
    the rest are untouched, so they may be contain dummy values sometimes!!!
=#
struct MetNet{T<:Real}
    N::Integer # number of fluxes
    M::Integer # number of metabolites
    S::AbstractMatrix{T} # Stoichiometric matrix M x N sparse (*)
    b::AbstractArray{T} # right hand side of equation  S ν = b (*)
    c::AbstractArray{Float64} # reaction index of biomass 
    lb::AbstractArray{T} # fluxes lower bound M elements vector (*)
    ub::AbstractArray{T} # fluxes upper bound M elements vector (*)
    genes::AbstractArray{String} # gene names 
    rxnGeneMat::AbstractMatrix{Float64} # 
    grRules::AbstractArray{String} # gene-reaction rule N elements vector of strings (and / or allowed)
    mets::AbstractArray{String} # metabolites short-name M elements (**)
    rxns::AbstractArray{String} # reactions short-name N elements (*)
    metNames::AbstractArray{String} # metabolites long-names M elements (**)
    metFormulas::AbstractArray{String} # metabolites formula M elements (**)
    rxnNames::AbstractArray{String} # reactions long-names N elements (**)
    rev::AbstractArray{Bool} # reversibility of reactions N elements
    subSystems::AbstractArray{String} # cellular component of fluxes N elements
end