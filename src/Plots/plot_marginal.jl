const ep_color = :red
const fba_color = :blue
const hr_color = :black

# Single out
function plot_marginal!(p, metnet::MetNet, out::Union{FBAout, EPout}, ider; 
    h = 1, legend = false, xlabel = "flx", ylabel = "pdf",
    color = :black, label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    av_ = av(metnet, out, ider);
    μ_ = μ(metnet, out, ider);
    σ_ = σ(metnet, out, ider);
    sσ_ = sqrt(σ_)
    tN = marginal(metnet, out, ider)
    local_max_ = μ_ <= lb_ ? lb_ : μ_ >= ub_ ? ub_ : μ_ # local max

    plot!(p, legend = legend, xlabel = xlabel, ylabel = ylabel)
    if sσ_ == 0.0 || !(-Inf < pdf(tN, local_max_) < Inf)
        # print("zero std, type ($(typeof(out)))!!!")
        plot!(p, [lb_ - ub_/10, ub_ + ub_/10], [0.0, 0.0]; 
            label = "", color = color, kwargs...)
        plot!(p, [av_, av_], [0.0, h]; 
            label = label, color = color, kwargs...)
    else
        # print("non zero std, type ($(typeof(out)))!!!")
        Plots.plot!(p, x -> pdf(tN, x), lb_ - ub_/10, ub_ + ub_/10, lw = 5; 
            label = label, color = color, kwargs...)
    end

end


function plot_marginal!(p, metnet::MetNet, out::HRout, ider; 
    h = :ignore, 
    color = :black, label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    hist = hists(metnet, out, ider)
    plot!(p, normalize(hist, mode = :pdf); 
        label = label, color = color, kwargs...)
end

# 3 outs
function plot_marginal!(p, metnet::MetNet, 
        fbaout::FBAout, epout::EPout, hrout::HRout, 
        ider::Union{AbstractString, Integer}; 
        legend = false, xlabel = "flx", ylabel = "pdf",
        title = rxns(metnet, ider), kwargs...)
    
    plot!(p, title = title, legend = legend, xlabel = xlabel, ylabel = ylabel)
    pdf_maxval_ = pdf_maxval(metnet, [hrout, epout], ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    plot_marginal!(p, metnet, hrout, ider; 
        color = hr_color, label = "HR", kwargs...)
    plot_marginal!(p, metnet, epout, ider; color = ep_color, 
        h = pdf_maxval_ * 1.1, lw = 5, label = "EP", kwargs...)
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
    ls = :dash,  color = fba_color,  lw = 5, label = "FBA", kwargs...)
end

function plot_marginal(metnet::MetNet, 
        fbaout::FBAout, epout::EPout, hrout::HRout, 
        ider::Union{AbstractString, Integer};
        title = rxns(metnet, ider), 
        xlabel = "flx", ylabel = "pdf", kwargs...)

        p = plot(title = title, xlabel = xlabel, ylabel = ylabel)
        return plot_marginal!(p, metnet, fbaout, epout, hrout, ider; kwargs...)
end

plot_marginal(metnet::MetNet, 
        fbaout::FBAout, epout::EPout, hrout::HRout, 
        iders; kwargs...) = 
    [plot_marginal(metnet, fbaout, epout, hrout, ider; kwargs...) for ider in iders]



# 2 outs
function plot_marginal!(p, metnet::MetNet, 
    fbaout::FBAout, epout::EPout,
    ider::IDER_TYPE; kwargs...)

    pdf_maxval_ = pdf_maxval(metnet, [epout], ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    plot_marginal!(p, metnet, epout, ider; h = pdf_maxval_ * 1.1, 
        color = ep_color, lw = 5, label = "EP", kwargs...)
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
        ls = :dash,  color = fba_color,  lw = 5, label = "FBA", kwargs...)
end

function plot_marginal!(p, metnet::MetNet, 
    fbaout::FBAout, hrout::HRout,
    ider::IDER_TYPE; kwargs...)

    pdf_maxval_ = pdf_maxval(metnet, [hrout], ider)
    plot_marginal!(p, metnet, hrout, ider; 
        color = hr_color, label = "HR", kwargs...)
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
        ls = :dash,  color = fba_color,  lw = 5, label = "FBA", kwargs...)
end

function plot_marginal!(p, metnet::MetNet, 
        epout::EPout, hrout::HRout,
        ider::IDER_TYPE; kwargs...)

    pdf_maxval_ = pdf_maxval(metnet, [epout], ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    plot_marginal!(p, metnet, hrout, ider; h = pdf_maxval_ * 1.1, 
        color = hr_color, label = "HR", kwargs...)
    plot_marginal!(p, metnet, epout, ider; 
        color = ep_color, lw = 5, label = "EP", kwargs...)
end

function plot_marginal(metnet::MetNet, out1, out2,
        ider::IDER_TYPE;
        title = rxns(metnet, ider), xlabel = "flx", ylabel = "pdf",
        kwargs...)

    p = plot(title = title, xlabel = xlabel, ylabel = ylabel)
    return plot_marginal!(p, metnet, out1, out2, ider; kwargs...)
end

plot_marginal(metnet::MetNet, out1, out2,
        iders::Vector; kwargs...) = [plot_marginal(metnet, out1, out2, ider; kwargs...) for ider in iders]


function plot_marginal!(p, boundle::ChstatBoundle, ξ::Real, β::Real, 
    ider::IDER_TYPE; kwargs...)

    outs = []
    try push!(outs, get_fbaout(boundle, ξ)) catch KeyError end
    try push!(outs, get_epout(boundle, ξ, β)) catch KeyError end
    try push!(outs, get_hrout(boundle, ξ, β)) catch KeyError end
    metnet = get_metnet(boundle, ξ)
    pdf_maxval_ = pdf_maxval(metnet, outs, ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    length(outs) == 1 && plot_marginal!(p, metnet, outs[1], ider; h = pdf_maxval_ * 1.1, kwargs...)
    length(outs) == 2 && plot_marginal!(p, metnet, outs[1], outs[2], ider; h = pdf_maxval_ * 1.1, kwargs...)
    length(outs) == 3 && plot_marginal!(p, metnet, outs[1], outs[2], outs[3], ider; h = pdf_maxval_ * 1.1, kwargs...)
    return p
end

plot_marginal(boundle::ChstatBoundle, ξ::Real, β::Real, 
    ider::IDER_TYPE; 
    title = rxns(get_metnet(boundle, ξ), ider),
    kwargs...) = 
        plot_marginal!(plot(), boundle, ξ, β, ider; title = title,  kwargs...)

plot_marginal(boundle::ChstatBoundle, ξ::Real, β::Real, iders::Vector; kwargs...) = 
    [plot_marginal(boundle, ξ, β, ider; kwargs...) for ider in iders]

function plot_marginal_legend(ξ, β; digits = 3, lw = 3, legend = :best, kwargs...)
    ξ = round(ξ, digits = digits)
    β = round(β, digits = digits)

    p = Plots.plot(framestyle = :none, legend = legend)
    Plots.plot!(p, [0,0],[0,0], color = ep_color, label = "EP beta: $β, xi: $ξ", lw = lw)
    Plots.plot!(p, [0,0],[0,0], color = fba_color, label = "FBA xi: $ξ", lw = lw)
    Plots.plot!(p, [0,0],[0,0], color = hr_color, label = "HR", lw = lw)
end
