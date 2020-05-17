struct EPMat{T<:AbstractFloat} <: AbstractEPMat
    KK::AbstractArray{T,2}
    KKPD::AbstractArray{T,2}
    invKKPD::Matrix{T}
    KY::Vector{T}
    v::Vector{T}
    lb::Vector{T}
    ub::Vector{T}
end

function EPMat(K::AbstractArray{T}, Y::Vector{T}, lb::Vector{T}, ub::Vector{T}, beta::T) where T <: Real
    M,N = size(K)
    KKPD = Matrix(beta * K' * K)
    if beta != Inf
        return EPMat(copy(KKPD), copy(KKPD), zeros(T,N,N), beta * K' * Y, zeros(T,N),lb,ub)
    else
        error("I really should not be here")
    end
end