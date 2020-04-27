function drop_hrsamples!(bound::ChstatBoundle, ξ::Real, β::Real)
    old_hrout = get_hrout(bound, ξ, β)
    if has_samples(old_hrout)
        new_hrout = HRout(hrsamples(old_hrout), true)
        add_hrout!(bound, ξ, β, new_hrout)
    end
end

drop_hrsamples!(bound::ChstatBoundle, ξs::Vector, β::Real) = 
    [drop_hrsamples!(bound, ξ, β) for ξ in ξs]