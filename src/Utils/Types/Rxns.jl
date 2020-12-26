struct Rxn{T<:Real}
    id::String
    S::Vector{T}
    mets::Vector
    c::T
    lb::T
    ub::T

    function Rxn{T}(id::String; 
            mets = [],
            S = T[],
            c = zero(T),
            lb = zero(T),
            ub = zero(T)
        ) where {T<:Real}
        
        length(S) != length(mets) && error("'S' and 'mets' must have the same length")
        new{T}(id, S, mets, c, lb, ub)
    end
end

Rxn(id::String; kwargs...) = Rxn{Float64}(id; kwargs...)