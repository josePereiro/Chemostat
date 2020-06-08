function balance_str(boundle::ChstatBoundle, ξ::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key; digits = 50) 

    model = get_data(boundle, ξ, metnet_data_key)
    out = get_data(boundle, ξ, data_key)
    balance_str(model, out, ider; digits = digits)
end

function balance_str(boundle::ChstatBoundle, ξ::Real, β::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key; digits = 50) 

    model = get_data(boundle, ξ, metnet_data_key)
    out = get_data(boundle, ξ, β, data_key)
    balance_str(model, out, ider; digits = digits)
end