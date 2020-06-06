struct EPMat{T<:AbstractFloat} <: AbstractEPMat
    KK::AbstractArray{T,2}
    KKPD::AbstractArray{T,2}
    invKKPD::Matrix{T}
    KY::Vector{T}
    v::Vector{T}
    lb::Vector{T}
    ub::Vector{T}
end

function EPMat(K::AbstractArray{T}, Y::Vector{T}, lb::Vector{T}, ub::Vector{T}, alpha::T) where T <: Real
    M,N = size(K)
    KKPD = Matrix(alpha * K' * K)
    if alpha != Inf
        return EPMat(copy(KKPD), copy(KKPD), zeros(T,N,N), alpha * K' * Y, zeros(T,N),lb,ub)
    else
        error("I really should not be here")
    end
end