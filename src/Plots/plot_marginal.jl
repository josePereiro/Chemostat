const ep_color = :red
const fba_color = :blue
const hr_color = :black

# Single out
function plot_marginal!(p, metnet::MetNet, out::Union{FBAout, EPout}, ider; h = 1,
    color = :black, label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    av_ = av(metnet, out, ider);
    μ_ = μ(metnet, out, ider);
    σ_ = σ(metnet, out, ider);
    sσ_ = sqrt(σ_)
    tN = marginal(metnet, out, ider)
    local_max_ = μ_ <= lb_ ? lb_ : μ_ >= ub_ ? ub_ : μ_ # local max

    if sσ_ == 0.0 || !(-Inf < pdf(tN, local_max_) < Inf)
        # print("zero std, type ($(typeof(out)))!!!")
        plot!(p, [lb_ - ub_/10, ub_ + ub_/10], [0.0, 0.0]; 
            label = "", color = color, kwargs...)
        plot!(p, [av_, av_], [0.0, h]; 
            label = label, color = color, kwargs...)
    else
        # print("non zero std, type ($(typeof(out)))!!!")
        Plots.plot!(p, x -> pdf(tN, x), lb_ - ub_/10, ub_ + ub_/10, lw = 10; 
            label = label, color = color, kwargs...)
    end

end


function plot_marginal!(p, metnet::MetNet, out::HRout, ider; 
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
        ider::Union{AbstractString, Integer})

    pdf_maxval_ = pdf_maxval(metnet, [hrout, epout], ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    plot_marginal!(p, metnet, hrout, ider; color = hr_color, label = "HR")
    plot_marginal!(p, metnet, epout, ider; color = ep_color, 
        h = pdf_maxval_ * 1.1, lw = 10, label = "EP")
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
    ls = :dash,  color = fba_color,  lw = 5, label = "FBA")
end

function plot_marginal(metnet::MetNet, 
        fbaout::FBAout, epout::EPout, hrout::HRout, 
        ider::Union{AbstractString, Integer};
        title = rxns(metnet, ider), xlabel = "flx", ylabel = "pdf")

        p = plot(title = title, xlabel = xlabel, ylabel = ylabel)
        return plot_marginal!(p, metnet, fbaout, epout, hrout, ider)
end

plot_marginal(metnet::MetNet, 
        fbaout::FBAout, epout::EPout, hrout::HRout, 
        iders) = [plot_marginal(metnet, fbaout, epout, hrout, ider) for ider in iders]



# 2 outs
function plot_marginal!(p, metnet::MetNet, 
    fbaout::EPout, epout::HRout,
    ider::Union{AbstractString, Integer})

    pdf_maxval_ = pdf_maxval(metnet, [epout], ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    plot_marginal!(p, metnet, epout, ider; h = pdf_maxval_ * 1.1, 
        color = ep_color, lw = 10, label = "EP")
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
        ls = :dash,  color = fba_color,  lw = 5, label = "FBA")
end

function plot_marginal!(p, metnet::MetNet, 
    fbaout::FBAout, hrout::HRout,
    ider::Union{AbstractString, Integer})

    pdf_maxval_ = pdf_maxval(metnet, [hrout], ider)
    plot_marginal!(p, metnet, hrout, ider; color = hr_color, label = "HR")
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
        ls = :dash,  color = fba_color,  lw = 5, label = "FBA")
end

function plot_marginal!(p, metnet::MetNet, 
        epout::EPout, hrout::HRout,
        ider::Union{AbstractString, Integer})

    pdf_maxval_ = pdf_maxval(metnet, [epout], ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    plot_marginal!(p, metnet, hrout, ider; h = pdf_maxval_ * 1.1, 
        color = hr_color, label = "HR")
    plot_marginal!(p, metnet, epout, ider; color = ep_color, lw = 10, label = "EP")
end

function plot_marginal(metnet::MetNet, out1, out2,
        ider::Union{AbstractString, Integer},
        title = rxns(metnet, ider), xlabel = "flx", ylabel = "pdf")

    p = plot(title = title, xlabel = xlabel, ylabel = ylabel)
    return plot_marginal!(p, metnet, out1, out2, ider)
end

plot_marginal(metnet::MetNet, out1, out2,
        iders) = [plot_marginal(metnet, out1, out2, ider) for ider in iders]