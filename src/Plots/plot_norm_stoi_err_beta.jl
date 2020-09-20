function plot_norm_stoi_err_beta!(p, bundle::ChstatBundle, 
    ξ::Real, βs::Vector, ider::Union{AbstractString, Integer}; 
    norm_fun::Function = x -> maximum(abs.(x)),
    label = ider, lw = 3, kwargs...)

ξ = parse_ξ(bundle, ξ)
ξstr = round(ξ, digits = 3)

plot!(title = "xi: $ξstr", xlabel = "beta", ylabel = "normalized stoi err |S.v - b|")
labeled = false
try
    errs = []
    isolated_met = false
    for β in βs
        err = av_stoi_err_ep(bundle, ξ, β, ider)
        
        metnet = get_metnet(bundle, ξ)
        rxns = met_rxns(metnet, ider)
        if isempty(rxns) 
            isolated_met = true
            break
        end
        err = err / norm_fun(av_ep(bundle, ξ, β, rxns))

        push!(errs, err)
    end
    if !isolated_met
        plot!(p, βs, errs; label = label, lw = lw, kwargs...)
        labeled = true
    end
catch KeyError end

try
    label = !labeled ? label : ""
    errs = []
    isolated_met = false
    for β in βs
        err = av_stoi_err_hr(bundle, ξ, β, ider)
        
        metnet = get_metnet(bundle, ξ)
        rxns = met_rxns(metnet, ider)
        if isempty(rxns) 
            isolated_met = true
            break
        end
        err = err / norm_fun(av_hr(bundle, ξ, β, rxns))

        push!(errs, err)
    end
    if !isolated_met
        plot!(p, βs, errs; label = label, ls = :dash, lw = lw, kwargs...)
    end
catch KeyError end

end


function plot_norm_stoi_err_beta!(p, bundle::ChstatBundle, 
    ξ::Real, βs::Vector, iders::Vector; lw = 3, kwargs...) 

    colors = distinguishable_colors(length(iders))
    for (ider, color) in zip(iders, colors)
        plot_norm_stoi_err_beta!(p, bundle, ξ, βs, ider; 
            color = color, lw = lw, kwargs...)
    end
end