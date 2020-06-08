function collect_data(fun::Function, boundle::ChstatBoundle, 
        ξ::Real, β::Real, data_keys::Vector)
    data = []
    for key in data_keys
        if haskey(boundle, ξ, key)
            d = boundle[ξ, key]
            fun(d) && push!(data, d)
        end
        if haskey(boundle, ξ, β, key)
            d = boundle[ξ, β, key]
            fun(d) && push!(data, d)
        end
    end
    return data
end

collect_data(boundle::ChstatBoundle, ξ::Real, β::Real, data_keys::Vector) = 
    collect_data((x) -> true, boundle, ξ::Real, β::Real, data_keys::Vector)