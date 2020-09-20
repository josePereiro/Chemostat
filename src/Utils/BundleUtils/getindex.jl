function Base.getindex(bundle::ChstatBundle, ξ::Real, β::Real)
    k_tuple = (parse_ξ(bundle, ξ), parse_β(bundle, β))
    return bundle.data[k_tuple]
end

function Base.getindex(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol)
    data = bundle[ξ, β]
    return haskey(data, data_key) ? data[data_key] : error("key '$data_key' not present, current keys in data dict $(collect(keys(data)))")
end

Base.getindex(bundle::ChstatBundle, ξs::Vector, β::Real, data_key::Symbol) =
    [bundle[ξ, β, data_key] for ξ in ξs]

Base.getindex(bundle::ChstatBundle, ξi, β::Real, data_key::Symbol) =
    [bundle[ξ, β, data_key] for ξ in bundle.ξs[ξi]]

Base.getindex(bundle::ChstatBundle, ξ::Real, βs::Vector, data_key::Symbol) =
    [bundle[ξ, β, data_key] for β in βs]

Base.getindex(bundle::ChstatBundle, ξ::Real, βi, data_key::Symbol) =
    [bundle[ξ, β, data_key] for β in bundle.βs[βi]]

function Base.getindex(bundle::ChstatBundle, ξ::Real)
    ξ = parse_ξ(bundle, ξ)
    return bundle.data[ξ]
end

function Base.getindex(bundle::ChstatBundle, ξ::Real, data_key::Symbol)
    data = bundle[ξ]
    return haskey(data, data_key) ? data[data_key] : error("key '$data_key' not present, current keys in data dict $(collect(keys(data)))")
end

Base.getindex(bundle::ChstatBundle, ξs::Vector, data_key::Symbol) = 
    [bundle[ξ, data_key] for ξ in ξs]

Base.getindex(bundle::ChstatBundle, ξi, data_key::Symbol) = 
    [bundle[ξ, data_key] for ξ in bundle.ξs[ξi]]

Base.getindex(bundle::ChstatBundle, data_key::Symbol) = bundle.data[data_key]