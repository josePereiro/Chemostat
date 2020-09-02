struct Rxn{T<:Real}
    id::String
    S::Vector{T}
    mets::Vector
    c::T
    lb::T
    ub::T

    function Rxn{T}(id::String; kwargs...) where {T<:Real}
        kwargs = Dict(kwargs)
        S = get(kwargs, :S, T[])
        mets = get(kwargs, :rxns, [])
        length(S) != length(mets) && error("'S' and 'mets' must have the same length")

        new{T}(id, S, mets,
            get(kwargs, :c, zero(T)),
            get(kwargs, :lb, zero(T)),
            get(kwargs, :ub, zero(T)),
        )
    end
end

Rxn(id::String; kwargs...) = Rxn{Float64}(id; kwargs...)