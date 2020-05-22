function Base.haskey(boundle::ChstatBoundle, ξ::Real, data_key::Symbol)
    try
        return Base.haskey(get_data(boundle, ξ), data_key)
    catch KeyError end
    return false
end

Base.haskey(boundle::ChstatBoundle, ξs::Vector, data_key::Symbol) =
     all([haskey(boundle, ξ, data_key) for ξ in ξs])

function Base.haskey(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol)
    try
        return Base.haskey(get_data(boundle, ξ, β), data_key)
    catch KeyError end
    return false
    
end

Base.haskey(boundle::ChstatBoundle, ξ, β, data_key::Symbol) =
    all([all([haskey(boundle, ξ_, β_, data_key) for β_ in β]) for ξ_ in ξ])