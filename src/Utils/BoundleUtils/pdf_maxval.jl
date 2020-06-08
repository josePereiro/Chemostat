function pdf_maxval(boundle::ChstatBoundle, 
            ξ::Real, β::Real, data_keys::Vector, 
            ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key)
        outs = collect_data(d -> d isa AbstractOut, boundle, ξ, β, data_keys)
        metnet = boundle[ξ, metnet_data_key]
        return pdf_maxval(metnet, outs, ider)
end