function Base.getindex(boundle::ChstatBoundle, ξ::Real, β::Real)
    k_tuple = (parse_ξ(boundle, ξ), parse_β(boundle, β))
    return boundle.data[k_tuple]
end

function Base.getindex(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol)
    data = boundle[ξ, β]
    return haskey(data, data_key) ? data[data_key] : error("key '$data_key' not present, current keys in data dict $(collect(keys(data)))")
end

Base.getindex(boundle::ChstatBoundle, ξs::Vector, β::Real, data_key::Symbol) =
    [boundle[ξ, β, data_key] for ξ in ξs]

Base.getindex(boundle::ChstatBoundle, ξi, β::Real, data_key::Symbol) =
    [boundle[ξ, β, data_key] for ξ in boundle.ξs[ξi]]

Base.getindex(boundle::ChstatBoundle, ξ::Real, βs::Vector, data_key::Symbol) =
    [boundle[ξ, β, data_key] for β in βs]

Base.getindex(boundle::ChstatBoundle, ξ::Real, βi, data_key::Symbol) =
    [boundle[ξ, β, data_key] for β in boundle.βs[βi]]

function Base.getindex(boundle::ChstatBoundle, ξ::Real)
    ξ = parse_ξ(boundle, ξ)
    return boundle.data[ξ]
end

function Base.getindex(boundle::ChstatBoundle, ξ::Real, data_key::Symbol)
    data = boundle[ξ]
    return haskey(data, data_key) ? data[data_key] : error("key '$data_key' not present, current keys in data dict $(collect(keys(data)))")
end

Base.getindex(boundle::ChstatBoundle, ξs::Vector, data_key::Symbol) = 
    [boundle[ξ, data_key] for ξ in ξs]

Base.getindex(boundle::ChstatBoundle, ξi, data_key::Symbol) = 
    [boundle[ξ, data_key] for ξ in boundle.ξs[ξi]]

Base.getindex(boundle::ChstatBoundle, data_key::Symbol) = boundle.data[data_key]