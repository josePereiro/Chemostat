norm_stoi_err(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol, 
            metnet_data_key::Symbol = metnet_data_key; normfun = norm1_stoi_err) = 
            normfun(bundle[ξ, metnet_data_key], bundle[ξ, β, data_key])

norm_stoi_err(bundle::ChstatBundle, ξ::Real, data_key::Symbol, 
        metnet_data_key::Symbol = metnet_data_key; normfun = norm1_stoi_err) = 
        normfun(bundle[ξ, metnet_data_key], bundle[ξ, data_key])

norm_stoi_err(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key; normfun = norm1_stoi_err) = 
        normfun(bundle[ξ, metnet_data_key], bundle[ξ, β, data_key],  ider)

norm_stoi_err(bundle::ChstatBundle, ξ::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key; normfun = norm1_stoi_err) = 
        normfun(bundle[ξ, metnet_data_key], bundle[ξ, data_key],  ider)