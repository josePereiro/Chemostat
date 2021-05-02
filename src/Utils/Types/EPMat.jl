struct EPMat{T<:AbstractFloat} <: AbstractEPMat
    αStS::SparseMatrixCSC{T,Int64} # A container of KK = βS'S (Its diag will change)
    αStS_diag::AbstractArray{T} # Store the original diag of KK = βS'S
    αStb::AbstractArray{T}
    Σ::Matrix{T} # ΣQ Enforced to be dense 
    v::AbstractArray{T}
    lb::AbstractArray{T}
    ub::AbstractArray{T}
end

function EPMat(S::AbstractArray{T}, b::AbstractArray{T}, 
        lb::AbstractArray{T}, ub::AbstractArray{T}, alpha::T) where T <: Real
    isinf(alpha) && error("For Inf alpha use 'EPMatT0'")
    M, N = size(S)
    αStS = sparse(alpha * S' * S)
    return EPMat(   
        #= αStS      =# αStS, 
        #= αStS_diag =# diag(αStS), 
        #= αStb      =# alpha * S' * b, 
        #= Σ         =# zeros(T,N,N), 
        #= v         =# zeros(T,N), 
        #= lb        =# lb, 
        #= ub        =# ub
    )
end