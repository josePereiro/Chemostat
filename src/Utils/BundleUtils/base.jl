function Base.haskey(bundle::ChstatBundle, ξ::Real, data_key::Symbol)
    try
        return Base.haskey(bundle[ξ], data_key)
    catch KeyError end
    return false
end

Base.haskey(bundle::ChstatBundle, ξs::Vector, data_key::Symbol) =
     all([haskey(bundle, ξ, data_key) for ξ in ξs])

function Base.haskey(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol)
    try
        return Base.haskey(bundle[ξ, β], data_key)
    catch KeyError end
    return false
    
end

Base.haskey(bundle::ChstatBundle, ξ, β, data_key::Symbol) =
    all([all([haskey(bundle, ξ_, β_, data_key) for β_ in β]) for ξ_ in ξ])