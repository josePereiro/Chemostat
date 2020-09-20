function drop_hrsamples!(bundle::ChstatBundle, ξ::Real, β::Real, hrkey = hrout_data_key)
    old_hrout = bundle[ξ, β, hrkey]
    if has_samples(old_hrout)
        new_hrout = HRout(hrsamples(old_hrout), true)
        bundle[ξ, β, hrkey] = new_hrout
    end
end

drop_hrsamples!(bundle::ChstatBundle, ξs::Vector, β::Real, hrkey = hrout_data_key) = 
    [drop_hrsamples!(bundle, ξ, β, hrkey) for ξ in ξs]