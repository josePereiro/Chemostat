function pdf_maxval(bundle::ChstatBundle, 
            ξ::Real, β::Real, data_keys::Vector, 
            ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key)
        outs = collect_data(d -> d isa AbstractOut, bundle, ξ, β, data_keys)
        metnet = bundle[ξ, metnet_data_key]
        return pdf_maxval(metnet, outs, ider)
end