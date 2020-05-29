stoi_err(bound::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol, 
            metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(get_data(bound, ξ, metnet_data_key), get_data(bound, ξ, β, data_key))


stoi_err(bound::ChstatBoundle, ξ::Real, data_key::Symbol, 
        metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(get_data(bound, ξ, metnet_data_key), get_data(bound, ξ, data_key))



stoi_err(bound::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(get_data(bound, ξ, metnet_data_key), get_data(bound, ξ, β, data_key),  ider)

stoi_err(bound::ChstatBoundle, ξ::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key) = 
    stoi_err(get_data(bound, ξ, metnet_data_key), get_data(bound, ξ, data_key),  ider)