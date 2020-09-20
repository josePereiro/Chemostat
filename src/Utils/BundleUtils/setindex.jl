function Base.setindex!(bundle::ChstatBundle, data, ξ::Real, β::Real, data_key::Symbol)
    _check_data_key(data_key, data)
    k_tuple = (ξ, β)
    if !haskey(bundle.data, k_tuple)
        bundle.data[k_tuple] = Dict()
    end
    bundle.data[k_tuple][data_key] = data

    !(ξ in bundle.ξs) && push!(bundle.ξs, ξ)
    !(β in bundle.βs) && push!(bundle.βs, β)
end
 
function Base.setindex!(bundle::ChstatBundle, data, ξ::Real, data_key::Symbol)
    _check_data_key(data_key, data)
    if !haskey(bundle.data, ξ)
        bundle.data[ξ] = Dict()
    end
    bundle.data[ξ][data_key] = data

    !(ξ in bundle.ξs) && push!(bundle.ξs, ξ)
end

function Base.setindex!(bundle::ChstatBundle, data, data_key::Symbol)
    _check_data_key(data_key, data)
    if !haskey(bundle.data, data_key)
        bundle.data[data_key] = Dict()
    end
    bundle.data[data_key] = data
end

# Warn about using non default keys for common data
function _check_data_key(data_key, data)
    _warn(default_key) = 
        data_key != default_key && @warn("Not using default data_key ($default_key)!!!")
#     data isa EPout ? _warn(epout_data_key) :
#     data isa FBAout ? _warn(fbaout_data_key) :
#     data isa HRout ? _warn(hrout_data_key) :
    data isa MetNet ? _warn(metnet_data_key) : nothing
end