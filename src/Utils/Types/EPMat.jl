struct EPMat{T<:AbstractFloat} <: AbstractEPMat
    KKc::AbstractArray{T,2} # A container of KK = βF'F (Its diag will change)
    dKKbk::AbstractVector{T} # Store the original diag of KK = βF'F
    invKKPD::Matrix{T}
    KY::Vector{T}
    v::Vector{T}
    lb::Vector{T}
    ub::Vector{T}
end

function EPMat(K::AbstractArray{T}, Y::Vector{T}, lb::Vector{T}, ub::Vector{T}, alpha::T) where T <: Real
    alpha == Inf && error("For Inf alpha use 'EPMatT0'")
    M,N = size(K)
    KKc = Matrix(alpha * K' * K)
    return EPMat(copy(KKc), diag(KKc), zeros(T,N,N), alpha * K' * Y, zeros(T,N), lb, ub)
end