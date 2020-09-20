function balance_str(bundle::ChstatBundle, ξ::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key; digits = 50) 

    model = bundle[ξ, metnet_data_key]
    out = bundle[ξ, data_key]
    balance_str(model, out, ider; digits = digits)
end

function balance_str(bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol, ider::IDER_TYPE,
        metnet_data_key::Symbol = metnet_data_key; digits = 50) 

    model = bundle[ξ, metnet_data_key]
    out = bundle[ξ, β, data_key]
    balance_str(model, out, ider; digits = digits)
end