const ep_color = :red
const fba_color = :blue
const hr_color = :black

function plot_marginal!(p, μ::Real, σ::Real, lb::Real, ub::Real;
        h = 1, lw = 5,
        color = :black,
        label = "", kwargs...)

    sσ = sqrt(σ)
    tN = Truncated(Normal(μ, sσ), lb, ub) 
    global_max_ = μ <= lb ? lb : μ >= ub ? ub : μ 
    margin_ = abs(ub - lb) * 0.1
    if sσ == 0.0 || !(-Inf < pdf(tN, global_max_) < Inf)
        Plots.plot!(p, [lb - margin_, ub + margin_], [0.0, 0.0]; 
            label = "")
        Plots.plot!(p, [av, av], [0.0, h]; 
            label = label, color = color, lw = lw, kwargs...)
    else
        Plots.Plots.plot!(p, x -> pdf(tN, x), lb - margin_, ub + margin_; 
            lw = lw, label = label, color = color, kwargs...)
    end
    return Plots.plot!(p; kwargs...)
end

# Single out
function plot_marginal!(p, metnet::MetNet, out::Union{FBAout, EPout}, ider; 
        color = out isa FBAout ? fba_color : ep_color,
        label = rxns(metnet, ider), kwargs...)
        
    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    μ_ = μ(metnet, out, ider);
    σ_ = σ(metnet, out, ider);
    plot_marginal!(p, μ_, σ_, lb_, ub_; color = color, label = label, kwargs...)
end


function plot_marginal!(p, metnet::MetNet, out::HRout, ider; 
        h = :ignored, color = hr_color, 
        label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)
    margin_ = abs(ub_ - lb_) * 0.1
    hist = hists(metnet, out, ider)
    Plots.plot!(p, normalize(hist, mode = :pdf); 
        xaxis = [lb_ - margin_, ub_ + margin_],
        label = label, color = color, kwargs...)
end

function plot_marginal!(p, metnet::MetNet, outs::Vector, ider::IDER_TYPE; kwargs...)
    pdf_maxval_ = pdf_maxval(metnet, outs, ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    for out in outs
        plot_marginal!(p, metnet, out, ider; h = pdf_maxval_ * 1.1, kwargs...)
    end
    return p
end

# Boundle
function plot_marginal!(p, boundle::ChstatBoundle, 
        ξ::Real, β::Real, data_keys::Vector,
        ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key; 
        kwargs...)

    metnet = boundle[ξ, metnet_data_key]
    outs = collect_data(d -> d isa AbstractOut, boundle, ξ, β, data_keys)
    plot_marginal!(p, metnet, outs, ider; kwargs...)
    return p
end