function get_data(boundle::ChstatBoundle, ξ::Real, β::Real)
    k_tuple = (parse_ξ(boundle, ξ), parse_β(boundle, β))
    return boundle.data[k_tuple]
end

function get_data(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol)
    data = get_data(boundle, ξ, β)
    return haskey(data, data_key) ? data[data_key] : error("key '$data_key' not present, currect keys in data dict $(collect(keys(data)))")
end

get_data(boundle::ChstatBoundle, ξs::Vector, β::Real, data_key::Symbol) =
    [get_data(boundle, ξ, β, data_key) for ξ in ξs]

get_data(boundle::ChstatBoundle, ξ::Real, βs::Vector, data_key::Symbol) =
    [get_data(boundle, ξ, β, data_key) for β in βs]

function get_data(boundle::ChstatBoundle, ξ::Real)
    ξ = parse_ξ(boundle, ξ)
    return boundle.data[ξ]
end

function get_data(boundle::ChstatBoundle, ξ::Real, data_key::Symbol)
    data = get_data(boundle, ξ)
    return haskey(data, data_key) ? data[data_key] : error("key '$data_key' not present, current keys in data dict $(collect(keys(data)))")
end

get_data(boundle::ChstatBoundle, ξs::Vector, data_key::Symbol) = 
    [get_data(boundle, ξ, data_key) for ξ in ξs]

get_data(boundle::ChstatBoundle, data_key::Symbol) = boundle.data[data_key]


