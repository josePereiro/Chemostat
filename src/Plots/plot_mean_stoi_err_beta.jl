function plot_mean_stoi_err_beta!(p, boundle::ChstatBoundle, 
    ξ::Real, βs::Vector; lw = 3, kwargs...) 

    ξ = parse_ξ(boundle, ξ)
    ξstr = round(ξ, digits = 3)
    plot!(title = "xi: $ξstr", xlabel = "beta", ylabel = "mean stoi err |S.v - b|")
    try
        errs = [mean(av_stoi_err_ep(boundle, ξ, β)) for β in βs]
        plot!(p, βs, errs; label = "EP", lw = lw, kwargs...)
    catch KeyError end
    try
        errs = [mean(av_stoi_err_hr(boundle, ξ, β)) for β in βs]
        plot!(p, βs, errs; label = "HR", ls = :dash, lw = lw, kwargs...)
    catch KeyError end

end