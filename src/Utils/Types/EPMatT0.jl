# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

struct EPMatT0{T<:AbstractFloat} <: AbstractEPMat
    Σy::Matrix{T}
    Σw::Matrix{T}
    G::Matrix{T}
    lb::Vector{T}
    ub::Vector{T}
    vy::Vector{T}
    vw::Vector{T}
    Y::Vector{T}
    idx::Vector{Int}
end

function EPMatT0(K::AbstractArray{T,2}, Y::Vector{T}, lb::Vector{T}, ub::Vector{T}) where T <: Real
    M,N = size(K)    
    M <= N || @warn("numeber of rows M=$M larger than number of cols N=$N")
    _,ci,EK,EY = echelonize(Matrix(K),Y)
    Mech,Nech = size(EK)
    length(EY) == Mech || error("vector size incompatible with matrix") 
    return EPMatT0(zeros(T,Mech,Mech), zeros(T,Nech-Mech,Nech-Mech), copy(EK[1:Mech,Mech+1:Nech]), lb[ci],ub[ci],zeros(T,Mech),zeros(T,Nech-Mech),EY,sortperm(ci)) # soerperm ci is the inverse permutation that sends me back to the original permutation
end