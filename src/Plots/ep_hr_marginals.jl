# Hit and Run marginal
function plot_marginal!(p, model::MetNet, hrsamples::AbstractArray, 
        ider::Union{Int, AbstractString}; kwargs...)
    kwargs = Dict(kwargs)
    ider = rxnindex(model, ider)
    samples = hrsamples[:,ider]
    av = round(mean(samples), digits = 2)
    va = round(var(samples), digits = 2)
    p = plot!(p, xlabel = "flx", ylabel = "pdf", 
        legend = get(kwargs, :legend, false), 
        title = get(kwargs, :title, "$(model.rxns[ider]) ~ (av: $av, va: $va)"))
    return Plots.histogram!(p, samples, normalize = :pdf, 
        color = get(kwargs, :color, :black), 
        label = get(kwargs, :label, "HR ($(model.rxns[ider]))"))
end
plot_marginal(model::MetNet, hrsamples::AbstractArray, ider::Union{Int, AbstractString}; kwargs...) =
    plot_marginal!(plot(), model, hrsamples, ider; kwargs...)

plot_marginals(model::MetNet, hrsamples::AbstractArray, iders; kwargs...) = 
    [plot_marginal(model, hrsamples, ider; kwargs...) for ider in iders]

# EPuot marginal
function plot_marginal!(p, model::MetNet, epout::EPout, ider::Union{Int, AbstractString}; kwargs...)
    ider = rxnindex(model, ider)
    tn = trunc_normal(model, epout, ider)
    av = round(mean(tn), digits = 2)
    va = round(var(tn), digits = 2)
    p = plot!(p, xlabel = "flx", ylabel = "pdf",
        legend = get(kwargs, :legend, false), 
        title = get(kwargs, :title, "$(model.rxns[ider]) ~ tN(av: $av, va: $va)"))
    
    lb = model.lb[ider]; ub = model.ub[ider]
    return plot!(p, x-> pdf(tn, x), lb - ub/10, ub + ub/10, 
        lw = get(kwargs, :lw, 10), 
        color = get(kwargs, :color, :red),
        label = get(kwargs, :label, "EP ($(model.rxns[ider]))"))
end
plot_marginal(model, epout::EPout, ider::Union{Int, AbstractString}; kwargs...) =
    plot_marginal!(plot(), model, epout, ider; kwargs...)

plot_marginals(model::MetNet, epout::EPout, iders; kwargs...) = 
    [plot_marginal(model, epout::EPout, ider; kwargs...) for ider in iders]

function plot_marginal(model::MetNet, hrsamples::AbstractArray, epout::EPout, 
        ider::Union{Int, AbstractString}; kwargs...)
    hist_ = plot_marginal(model, hrsamples, ider)
    return plot_marginal!(hist_, model, epout, ider; kwargs...)
end
plot_marginal(model::MetNet, epout::EPout, hrsamples::AbstractArray, 
            ider::Union{Int, AbstractString}; kwargs...) = 
        plot_marginal(model, hrsamples, epout, ider; kwargs...)

plot_marginals(model::MetNet, hrsamples::AbstractArray, epout::EPout, iders; kwargs...) = 
    [plot_marginal(model, hrsamples, epout, ider; kwargs...) for ider in iders]

plot_marginals(model::MetNet, epout::EPout, hrsamples::AbstractArray, iders; kwargs...) = 
    plot_marginals(model, hrsamples, epout, iders; kwargs...)