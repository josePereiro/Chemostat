function find_closest_beta(bundle::ChstatBundle, ξ::Real, target::Real, ider::IDER_TYPE)
    βs = bundle.βs
    exp_β = βs |> first
    for β in βs
        ep_μ = av(bundle, ξ, β, :ep, ider)
        last_ep_μ = av(bundle, ξ, exp_β, :ep, ider)
        if abs(ep_μ - target) < abs(last_ep_μ - target)
            exp_β = β
        end
    end
    return exp_β
end
