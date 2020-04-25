function plot_mean_stoi_err_xi!(p, boundle::ChstatBoundle, 
    ξs::Vector, β::Real; lw = 3, kwargs...) 

    β = parse_β(boundle, β)
    βstr = round(β, digits = 3)

    plot!(title = "beta: $βstr", xlabel = "xi", ylabel = "mean stoi err |S.v - b|")
    try
        errs = [mean(av_stoi_err_ep(boundle, ξ, β)) for ξ in ξs]
        plot!(p, ξs, errs; label = "EP", lw = lw, kwargs...)
    catch KeyError end
    try
        errs = [mean(av_stoi_err_hr(boundle, ξ, β)) for ξ in ξs]
        plot!(p, ξs, errs; label = "HR", ls = :dash, lw = lw, kwargs...)
    catch KeyError end

end