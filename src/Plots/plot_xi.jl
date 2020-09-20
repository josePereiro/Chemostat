function plot_ξs!(p, ξs::Vector, models::Vector, outs::Vector, 
        ider::Union{AbstractString, Integer}; 
        xlabel = "xi", ylabel = "flx",
        title = ider, lw = 3, color = :black, 
        legend = false, kwargs...)

    avs = av(models, outs, ider)
    vas = va(models, outs, ider)
    plot!(p, title = title, xlabel = xlabel, ylabel = ylabel, legend = legend)
    plot!(p, ξs, avs; lw = lw, color = color, kwargs...)
    plot!(p, ξs, sqrt.(vas); lw = lw, color = color, ls = :dash, kwargs...)
end

plot_ξs!(p, ξs::Vector, models::Vector, outs::Vector, 
    iders::Vector; lw = 3, color = :black, kwargs...) = 
        [plot_ξs!(p, ξs, models, outs, ider; lw = lw, color = color, kwargs...) 
            for ider in iders]

function plot_ξs!(p, ξs::Vector, β::Real, bundle::ChstatBundle, 
    ider::Union{AbstractString, Integer}; 
    lw = 3, kwargs...)

    try 
        plot_ξs!(p, ξs, get_metnet(bundle, ξs),
                get_epout(bundle, ξs, β), 
                ider; lw = lw, color = ep_color, label = "EP", kwargs...)
    catch KeyError end
    try 
        plot_ξs!(p, ξs, get_metnet(bundle, ξs),
                get_hrout(bundle, ξs, β), 
                ider; lw = lw, color = hr_color, label = "HR", kwargs...)
    catch KeyError end
    try 
        plot_ξs!(p, ξs, get_metnet(bundle, ξs),
                get_fbaout(bundle, ξs), 
                ider; lw = lw, color = fba_color, label = "FBA", kwargs...)
    catch KeyError end
end

plot_ξs(ξs::Vector, β::Real, bundle::ChstatBundle, 
    ider::Union{AbstractString, Integer}; 
    lw = 3, kwargs...) = 
        plot_ξs!(plot(), ξs, β, bundle, ider; lw = lw, kwargs...)

plot_ξs(ξs::Vector, β::Real, bundle::ChstatBundle, 
    iders::Vector; lw = 3, kwargs...) = 
        [plot_ξs(ξs, β, bundle, ider; lw = lw, kwargs...)
            for ider in iders]

function plot_ξs_legend(β; digits = 3, lw = 3, legend = :best, kwargs...)
    β = round(β, digits = digits)
    # legend
    p = plot(;framestyle = :none, legend = legend, kwargs...)
    plot!(p, [0,0],[0,0], color = fba_color, label = "ave FBA", lw = 3)
    plot!(p, [0,0],[0,0], color = ep_color, label = "ave EP beta: $β", lw = 3)
    plot!(p, [0,0],[0,0],  color = hr_color, label = "ave HR", lw = 3)
    plot!(p, [0,0],[0,0], color = fba_color, label = "std FBA", lw = 3, ls = :dash)
    plot!(p, [0,0],[0,0], color = ep_color, label = "std EP beta: $β", lw = 3, ls = :dash)
    plot!(p, [0,0],[0,0], color = hr_color, label = "std HR", lw = 3, ls = :dash)
end