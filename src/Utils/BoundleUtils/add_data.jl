function add_data!(boundle::ChstatBoundle, ξ::Real, β::Real, key::AbstractString, data)
    k_tuple = (ξ, β)
    if !haskey(boundle.data, k_tuple)
        boundle.data[k_tuple] = Dict()
    end
    boundle.data[k_tuple][key] = data

    !(ξ in boundle.ξs) && push!(boundle.ξs, ξ)
    !(β in boundle.βs) && push!(boundle.βs, β)
end

# add_epout!(boundle::ChstatBoundle, ξ::Real, β::Real, epout::EPout) = 
#     add_data!(boundle, ξ, β, epout_key, epout)

# add_epmat!(boundle::ChstatBoundle, ξ::Real, β::Real, epmat::AbstractEPMat) = 
#     add_data!(boundle, ξ, β, epmat_key, epmat)

# add_hrout!(boundle::ChstatBoundle, ξ::Real, β::Real, hrout::HRout) = 
#     add_data!(boundle, ξ, β, hrout_key, hrout)
 
function add_data!(boundle::ChstatBoundle, ξ::Real, key::AbstractString, data)
    if !haskey(boundle.data, ξ)
        boundle.data[ξ] = Dict()
    end
    boundle.data[ξ][key] = data

    !(ξ in boundle.ξs) && push!(boundle.ξs, ξ)
end

# add_fbaout!(boundle::ChstatBoundle, ξ::Real, fbaout::FBAout) = 
#     add_data!(boundle, ξ, fbaout_key, fbaout)

# add_metnet!(boundle::ChstatBoundle, ξ::Real, metnet::MetNet) = 
#     add_data!(boundle, ξ, metnet_key, metnet)

function add_data!(boundle::ChstatBoundle, key::AbstractString, data)
    if !haskey(boundle.data, key)
        boundle.data[key] = Dict()
    end
    boundle.data[key] = data
end