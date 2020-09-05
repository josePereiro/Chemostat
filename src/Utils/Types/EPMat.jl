struct EPMat{T<:AbstractFloat} <: AbstractEPMat
    KKc::SparseMatrixCSC{T,Int64} # A container of KK = βF'F (Its diag will change)
    dKKbk::AbstractArray{T} # Store the original diag of KK = βF'F
    invKKPD::Matrix{T} # Enforced to be dense 
    KY::AbstractArray{T}
    v::AbstractArray{T}
    lb::AbstractArray{T}
    ub::AbstractArray{T}
end

function EPMat(K::AbstractArray{T}, Y::AbstractArray{T}, lb::AbstractArray{T}, ub::AbstractArray{T}, alpha::T) where T <: Real
    alpha == Inf && error("For Inf alpha use 'EPMatT0'")
    M,N = size(K)
    Kt = K'
    KKc = sparse(alpha * Kt * K)
    return EPMat(KKc, diag(KKc), zeros(T,N,N), alpha * Kt * Y, zeros(T,N), lb, ub)
end