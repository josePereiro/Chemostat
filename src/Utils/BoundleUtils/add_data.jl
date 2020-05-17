function add_data!(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol, data)
    check_data_key(data_key, data)
    k_tuple = (ξ, β)
    if !haskey(boundle.data, k_tuple)
        boundle.data[k_tuple] = Dict()
    end
    boundle.data[k_tuple][data_key] = data

    !(ξ in boundle.ξs) && push!(boundle.ξs, ξ)
    !(β in boundle.βs) && push!(boundle.βs, β)
end
 
function add_data!(boundle::ChstatBoundle, ξ::Real, data_key::Symbol, data)
    check_data_key(data_key, data)
    if !haskey(boundle.data, ξ)
        boundle.data[ξ] = Dict()
    end
    boundle.data[ξ][data_key] = data

    !(ξ in boundle.ξs) && push!(boundle.ξs, ξ)
end

function add_data!(boundle::ChstatBoundle, data_key::Symbol, data)
    check_data_key(data_key, data)
    if !haskey(boundle.data, data_key)
        boundle.data[data_key] = Dict()
    end
    boundle.data[data_key] = data
end

# Warn about using non default keys for common data
function check_data_key(data_key, data)
    _warn(default_key) = 
        data_key != default_key && @warn("Not using default data_key ($default_key)!!!")
    data isa EPout ? _warn(epout_data_key) :
    data isa FBAout ? _warn(fbaout_data_key) :
    data isa HRout ? _warn(hrout_data_key) :
    data isa MetNet ? _warn(metnet_data_key) : nothing
end