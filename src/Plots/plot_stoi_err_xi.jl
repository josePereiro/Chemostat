function plot_stoi_err_xi!(p, bundle::ChstatBundle, 
        ξs::Vector, β::Real, ider::Union{AbstractString, Integer}; 
        label = ider, lw = 3, kwargs...)

    β = parse_β(bundle, β)
    βstr = round(β, digits = 3)

    plot!(title = "beta: $βstr", xlabel = "xi", ylabel = "stoi err |S.v - b|")
    labeled = false
    try
        errs = [av_stoi_err_ep(bundle, ξ, β, ider) for ξ in ξs]
        plot!(p, ξs, errs; label = label, lw = lw, kwargs...)
        labeled = true
    catch KeyError end
    try
        label = !labeled ? label : ""
        errs = [av_stoi_err_hr(bundle, ξ, β, ider) for ξ in ξs]
        plot!(p, ξs, errs; label = label, ls = :dash, lw = lw, kwargs...)
    catch KeyError end
end

function plot_stoi_err_xi!(p, bundle::ChstatBundle, 
        ξs::Vector, β::Real, iders; lw = 3, kwargs...) 
    
    colors = distinguishable_colors(length(iders))
    for (ider, color) in zip(iders, colors)
        plot_stoi_err_xi!(p, bundle, ξs, β, ider; 
            color = color, lw = lw, kwargs...)
    end
end