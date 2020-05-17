function drop_hrsamples!(bound::ChstatBoundle, ξ::Real, β::Real, hrkey = hrout_data_key)
    old_hrout = get_data(bound, ξ, β, hrkey)
    if has_samples(old_hrout)
        new_hrout = HRout(hrsamples(old_hrout), true)
        add_data!(bound, ξ, β, hrkey, new_hrout)
    end
end

drop_hrsamples!(bound::ChstatBoundle, ξs::Vector, β::Real, hrkey = hrout_data_key) = 
    [drop_hrsamples!(bound, ξ, β, hrkey) for ξ in ξs]