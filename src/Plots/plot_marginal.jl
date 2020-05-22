const ep_color = :red
const fba_color = :blue
const hr_color = :black

# Single out
function plot_marginal!(p, metnet::MetNet, out::Union{FBAout, EPout}, ider; 
        h = 1, 
        color = out isa FBAout ? fba_color : ep_color,
        label = rxns(metnet, ider), kwargs...)
        
    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    av_ = av(metnet, out, ider);
    μ_ = μ(metnet, out, ider);
    σ_ = σ(metnet, out, ider);
    sσ_ = sqrt(σ_)
    tN = marginal(metnet, out, ider)
    local_max_ = μ_ <= lb_ ? lb_ : μ_ >= ub_ ? ub_ : μ_ # local max

    if sσ_ == 0.0 || !(-Inf < pdf(tN, local_max_) < Inf)
        plot!(p, [lb_ - ub_/10, ub_ + ub_/10], [0.0, 0.0]; 
            label = "", color = color, lw = 5, kwargs...)
        plot!(p, [av_, av_], [0.0, h]; 
            label = label, color = color, lw = 5, kwargs...)
    else
        Plots.plot!(p, x -> pdf(tN, x), lb_ - ub_/10, ub_ + ub_/10; lw = 5, 
            label = label, color = color, kwargs...)
    end

end


function plot_marginal!(p, metnet::MetNet, out::HRout, ider; 
        h = :ignore, color = hr_color, 
        label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    hist = hists(metnet, out, ider)
    plot!(p, normalize(hist, mode = :pdf); 
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

# 
function plot_marginal!(p, boundle::ChstatBoundle, 
        ξ::Real, β::Real, data_keys::Vector,
        ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key; 
        kwargs...)

    outs = []
    metnet = get_data(boundle, ξ, metnet_data_key)
    for key in data_keys
        if haskey(boundle, ξ, key)
            data = get_data(boundle, ξ, key)
            data isa AbstractOut && push!(outs, data)
        end
        if haskey(boundle, ξ, β, key)
            data = get_data(boundle, ξ, β, key)
            data isa AbstractOut && push!(outs, data)
        end
    end
    plot_marginal!(p, metnet, outs, ider; kwargs...)
    return p
end