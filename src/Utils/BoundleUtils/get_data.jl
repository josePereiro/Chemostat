function get_data(boundle::ChstatBoundle, ξ::Real, β::Real)
    k_tuple = (parse_ξ(boundle, ξ), parse_β(boundle, β))
    return boundle.data[k_tuple]
end

function get_data(boundle::ChstatBoundle, ξ::Real, β::Real, key::AbstractString)
    k_tuple = (parse_ξ(boundle, ξ), parse_β(boundle, β))
    return boundle.data[k_tuple][key]
end

get_data(boundle::ChstatBoundle, ξs::Vector, β::Real, key::AbstractString) =
    [get_data(boundle, ξ, β, key) for ξ in ξs]

get_data(boundle::ChstatBoundle, ξ::Real, βs::Vector, key::AbstractString) =
    [get_data(boundle, ξ, β, key) for β in βs]

get_data(boundle::ChstatBoundle, ξ::Real, β::Real, keys::Vector) =
    [get_data(boundle, ξ, β, key) for key in keys]

get_epout(boundle::ChstatBoundle, ξ, β) = get_data(boundle, ξ, β, epout_key)
get_epmat(boundle::ChstatBoundle, ξ, β) = get_data(boundle, ξ, β, epmat_key)
get_hrout(boundle::ChstatBoundle, ξ, β) = get_data(boundle, ξ, β, hrout_key)


function get_data(boundle::ChstatBoundle, ξ::Real)
    ξ = parse_ξ(boundle, ξ)
    return boundle.data[ξ]
end

function get_data(boundle::ChstatBoundle, ξ::Real, key::AbstractString)
    ξ = parse_ξ(boundle, ξ)
    return boundle.data[ξ][key]
end

get_data(boundle::ChstatBoundle, ξs::Vector, key::AbstractString) = 
    [get_data(boundle, ξ, key) for ξ in ξs]

get_data(boundle::ChstatBoundle, ξ::Real, keys::Vector) = 
    [get_data(boundle, ξ, key) for key in keys]

get_fbaout(boundle::ChstatBoundle, ξ) = get_data(boundle, ξ, fbaout_key)
get_metnet(boundle::ChstatBoundle, ξ) = get_data(boundle, ξ, metnet_key)

get_data(boundle::ChstatBoundle, key::AbstractString) = boundle.data[key]
