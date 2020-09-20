function plot_stoi_err_beta!(p, bundle::ChstatBundle, 
        ξ::Real, βs::Vector, ider::Union{AbstractString, Integer}; 
        label = ider, lw = 3, kwargs...)
    
    ξ = parse_ξ(bundle, ξ)
    ξstr = round(ξ, digits = 3)

    plot!(title = "xi: $ξstr", xlabel = "beta", ylabel = "stoi err |S.v - b|")
    labeled = false
    try
        errs = [av_stoi_err_ep(bundle, ξ, β, ider) for β in βs]
        plot!(p, βs, errs; label = label, lw = lw, kwargs...)
        labeled = true
    catch KeyError end
    try
        label = !labeled ? label : ""
        errs = [av_stoi_err_hr(bundle, ξ, β, ider) for β in βs]
        plot!(p, βs, errs; label = label, ls = :dash, lw = lw, kwargs...)
    catch KeyError end
end
function plot_stoi_err_beta!(p, bundle::ChstatBundle, 
        ξ::Real, βs::Vector, iders::Vector; lw = 3, kwargs...) 

    colors = distinguishable_colors(length(iders))
    for (ider, color) in zip(iders, colors)
        plot_stoi_err_beta!(p, bundle, ξ, βs, ider; 
            color = color, lw = lw, kwargs...)
    end
end