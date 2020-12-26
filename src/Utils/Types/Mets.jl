struct Met{T<:Real}
    id::String
    S::Vector{T}
    rxns::Vector
    b::T

    function Met{T}(id::String; 
            S = T[],
            rxns = [],
            b = zero(T)
        ) where {T<:Real}
        
        length(S) != length(rxns) && error("'S' and 'rxns' must have the same length")
        return new{T}(id, S, rxns, b)
    end
end

Met(id::String; kwargs...) = Met{Float64}(id; kwargs...)