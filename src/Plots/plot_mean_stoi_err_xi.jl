function plot_mean_stoi_err_xi!(p, bundle::ChstatBundle, 
    ξs::Vector, β::Real; lw = 3, kwargs...) 

    β = parse_β(bundle, β)
    βstr = round(β, digits = 3)

    plot!(title = "beta: $βstr", xlabel = "xi", ylabel = "mean stoi err |S.v - b|")
    try
        errs = [mean(av_stoi_err_ep(bundle, ξ, β)) for ξ in ξs]
        plot!(p, ξs, errs; label = "EP", lw = lw, kwargs...)
    catch KeyError end
    try
        errs = [mean(av_stoi_err_hr(bundle, ξ, β)) for ξ in ξs]
        plot!(p, ξs, errs; label = "HR", ls = :dash, lw = lw, kwargs...)
    catch KeyError end

end