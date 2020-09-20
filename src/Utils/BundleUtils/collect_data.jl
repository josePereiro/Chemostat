function collect_data(fun::Function, bundle::ChstatBundle, 
        ξ::Real, β::Real, data_keys::Vector)
    data = []
    for key in data_keys
        if haskey(bundle, ξ, key)
            d = bundle[ξ, key]
            fun(d) && push!(data, d)
        end
        if haskey(bundle, ξ, β, key)
            d = bundle[ξ, β, key]
            fun(d) && push!(data, d)
        end
    end
    return data
end

collect_data(bundle::ChstatBundle, ξ::Real, β::Real, data_keys::Vector) = 
    collect_data((x) -> true, bundle, ξ::Real, β::Real, data_keys::Vector)