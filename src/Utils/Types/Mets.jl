struct Met{T<:Real}
    id::String
    S::Vector{T}
    rxns::Vector
    b::T

    function Met{T}(id::String; kwargs...) where {T<:Real}
        kwargs = Dict(kwargs)
        
        S = get(kwargs, :S, T[])
        rxns = get(kwargs, :rxns, [])
        length(S) != length(rxns) && error("'S' and 'rxns' must have the same length")

        return new{T}(id, S, rxns, get(kwargs, :b, zero(T)))
    end
end

Met(id::String; kwargs...) = Met{Float64}(id; kwargs...)