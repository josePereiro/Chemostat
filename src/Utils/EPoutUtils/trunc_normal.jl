function trunc_normal(model::MetNet, epout::EPout, ider::Union{Int, AbstractString})
    ider = rxnindex(model, ider)
    return Truncated(Normal(epout.μ[ider], sqrt(epout.σ[ider])), model.lb[ider], model.ub[ider]) 
end