# Boundle
av(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol) = 
    av(get_data(boundle, ξ, β, data_key))

function av(boundle::ChstatBoundle, ξ, β, data_key::Symbol, ider,
        metnet_data_key::Symbol = metnet_data_key)

    metnet = get_data(boundle, ξ, metnet_data_key)
    data = get_data(boundle, ξ, β, data_key)
    return av(metnet, data, ider)
end

va(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol) = 
    va(get_data(boundle, ξ, β, data_key))

function va(boundle::ChstatBoundle, ξ, β, data_key::Symbol, ider,
        metnet_data_key::Symbol = metnet_data_key)

    metnet = get_data(boundle, ξ, metnet_data_key)
    data = get_data(boundle, ξ, β, data_key)
    return va(metnet, data, ider)
end

μ(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol) = 
    μ(get_data(boundle, ξ, β, data_key))

function μ(boundle::ChstatBoundle, ξ, β, data_key::Symbol, ider,
        metnet_data_key::Symbol = metnet_data_key)

    metnet = get_data(boundle, ξ, metnet_data_key)
    data = get_data(boundle, ξ, β, data_key)
    return μ(metnet, data, ider)
end
σ(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol) = 
σ(get_data(boundle, ξ, β, data_key))

function σ(boundle::ChstatBoundle, ξ, β, data_key::Symbol, ider,
        metnet_data_key::Symbol = metnet_data_key)

    metnet = get_data(boundle, ξ, metnet_data_key)
    data = get_data(boundle, ξ, β, data_key)
    return σ(metnet, data, ider)
end