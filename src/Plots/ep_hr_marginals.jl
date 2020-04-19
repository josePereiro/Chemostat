function plot_marginal(epout, hrsamples, model, ider::Union{Int, AbstractString})
    ider = Utils.rxnindex(model, ider)
    p_ = plot(xlabel = "flx", ylabel = "pdf", legend = false)
    histogram!(p_, hrsamples[:,ider], normalize = :pdf)
    lb = model.lb[ider]
    ub = model.ub[ider]
    dist_ = Truncated(Normal(epout.μ[ider], sqrt(epout.σ[ider])), lb, ub)
    mean_str = string(round(mean(dist_), digits = 2))
    var_str = string(round(var(dist_), digits = 2))
    return plot!(p_, x-> pdf(dist_, x), lb - ub/10, ub + ub/10, lw = 10, 
        title = "$(model.rxns[ider]) ~ N($mean_str, $var_str)")
end

plot_marginals(epout, hrsamples, model, iders) = 
    [plot_marginal(epout, hrsamples, model, ider) for ider in iders]