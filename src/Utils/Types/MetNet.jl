# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

#=
    The * mark the relevant fields for EP algorithm, 
    and the ** the relevante for other methods in this package,
    the rest are untouched, so they may be contain dummy values sometimes!!!
=#
struct MetNet
    N::Int # number of fluxes
    M::Int # number of metabolites
    S::SparseMatrixCSC{Float64,Int} # Stoichiometric matrix M x N sparse (*)
    b::Array{Float64,1} # right hand side of equation  S ν = b (*)
    c::Array{Float64,1} # reaction index of biomass 
    lb::Array{Float64,1} # fluxes lower bound M elements vector (*)
    ub::Array{Float64,1} # fluxes upper bound M elements vector (*)
    genes::Array{String,1} # gene names 
    rxnGeneMat::SparseMatrixCSC{Float64,Int} # 
    grRules::Array{String,1} # gene-reaction rule N elements vector of strings (and / or allowed)
    mets::Array{String,1} # metabolites short-name M elements (**)
    rxns::Array{String,1} # reactions short-name N elements (*)
    metNames::Array{String,1} # metabolites long-names M elements (**)
    metFormulas::Array{String,1} # metabolites formula M elements (**)
    rxnNames::Array{String,1} # reactions long-names N elements (**)
    rev::Array{Bool,1} # reversibility of reactions N elements
    subSystems::Array{String,1} # cellular component of fluxes N elements
end;