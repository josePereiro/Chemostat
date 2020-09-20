stoi_err(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol, 
            metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(bundle[ξ, metnet_data_key], bundle[ξ, β, data_key])


stoi_err(bundle::ChstatBundle, ξ::Real, data_key::Symbol, 
        metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(bundle[ξ, metnet_data_key], bundle[ξ, data_key])



stoi_err(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(bundle[ξ, metnet_data_key], bundle[ξ, β, data_key],  ider)

stoi_err(bundle::ChstatBundle, ξ::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(bundle[ξ, metnet_data_key], bundle[ξ, data_key],  ider)