const ep_color = :red
const fba_color = :blue
const hr_color = :black


function plot_marginal!(p, metnet::MetNet, out::Union{FBAout, EPout}, ider; h = 1,
    color = :black, label = rxns(metnet, ider), kwargs...)
    kwargs = Dict(kwargs)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    dist = marginal(metnet, out, ider)
    if var(dist) < (ub_ - lb_)/50
        plot!(p, [lb_ - ub_/10, ub_ + ub_/10], [0.0, 0.0]; 
            label = "", color = color, kwargs...)
        plot!(p, [mean(dist), mean(dist)], [0.0, h]; 
            label = label, color = color, kwargs...)
    else
        plot!(p, x-> pdf(dist, x), lb_ - ub_/10, ub_ + ub_/10;
                    label = label, color = color, kwargs...)
    end

end

function plot_marginal!(p, metnet::MetNet, out::HRout, ider; 
    color = :black, label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    hist = hists(metnet, out, ider)
    plot!(p,normalize(hist, mode = :pdf); 
        label = label, color = color, kwargs...)
end

function plot_marginal!(p, metnet::MetNet, 
        hrout::HRout, epout::EPout, fbaout::FBAout, ider)

    pdf_maxval_ = pdf_maxval(metnet, [hrout, epout], ider)
    plot_marginal!(p, metnet, hrout, ider; color = hr_color, label = "HR")
    plot_marginal!(p, metnet, epout, ider; color = ep_color, lw = 10, label = "EP")
    plot_marginal!(p, metnet, fbaout, ider; h = pdf_maxval_ * 1.1,
    ls = :dash,  color = fba_color,  lw = 5, label = "FBA")
end

function plot_marginal(metnet::MetNet, hrout::HRout, epout::EPout, fbaout::FBAout, 
        ider::Union{AbstractString, Integer})
    p = Plots.plot(title = "$(rxns(metnet, ider))", xlabel = "flx", ylabel = "pdf")
    plot_marginal!(p, metnet, hrout, epout, fbaout, ider)
    return p
end

plot_marginal(metnet::MetNet, hrout::HRout, epout::EPout, fbaout::FBAout, iders) = 
    [plot_marginal(metnet, hrout, epout, fbaout, ider) for ider in iders]

plot_marginal(metnet::MetNet, hrout::HRout, fbaout::FBAout, epout::EPout, iders) = 
    [plot_marginal(metnet, hrout, epout, fbaout, ider) for ider in iders]

plot_marginal(metnet::MetNet, fbaout::FBAout, epout::EPout, hrout::HRout, iders) = 
    [plot_marginal(metnet, hrout, epout, fbaout, ider) for ider in iders]

plot_marginal(metnet::MetNet, fbaout::FBAout, hrout::HRout, epout::EPout, iders) = 
    [plot_marginal(metnet, hrout, epout, fbaout, ider) for ider in iders]

plot_marginal(metnet::MetNet, epout::EPout, fbaout::FBAout, hrout::HRout, iders) = 
    [plot_marginal(metnet, hrout, epout, fbaout, ider) for ider in iders]

plot_marginal(metnet::MetNet, epout::EPout, hrout::HRout, fbaout::FBAout, iders) = 
    [plot_marginal(metnet, hrout, epout, fbaout, ider) for ider in iders]