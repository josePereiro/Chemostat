norm_abs_stoi_err(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol, 
            metnet_data_key::Symbol = metnet_data_key) = 
    norm_abs_stoi_err(boundle[ξ, metnet_data_key], boundle[ξ, β, data_key])


norm_abs_stoi_err(boundle::ChstatBoundle, ξ::Real, data_key::Symbol, 
        metnet_data_key::Symbol = metnet_data_key) = 
    norm_abs_stoi_err(boundle[ξ, metnet_data_key], boundle[ξ, data_key])



norm_abs_stoi_err(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key) = 
    norm_abs_stoi_err(boundle[ξ, metnet_data_key], boundle[ξ, β, data_key],  ider)

norm_abs_stoi_err(boundle::ChstatBoundle, ξ::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key) = 
    norm_abs_stoi_err(boundle[ξ, metnet_data_key], boundle[ξ, data_key],  ider)