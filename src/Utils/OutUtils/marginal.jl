function marginal(model::MetNet, out::Union{FBAout, EPout}, ider)
    ider = rxnindex(model, ider)
    μ_ = μ(model, out, ider)
    σ_ = σ(model, out, ider)
    return Truncated(Normal(μ_, sqrt(σ_)), model.lb[ider], model.ub[ider]) 
end