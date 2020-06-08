function drop_hrsamples!(boundle::ChstatBoundle, ξ::Real, β::Real, hrkey = hrout_data_key)
    old_hrout = boundle[ξ, β, hrkey]
    if has_samples(old_hrout)
        new_hrout = HRout(hrsamples(old_hrout), true)
        boundle[ξ, β, hrkey] = new_hrout
    end
end

drop_hrsamples!(boundle::ChstatBoundle, ξs::Vector, β::Real, hrkey = hrout_data_key) = 
    [drop_hrsamples!(boundle, ξ, β, hrkey) for ξ in ξs]