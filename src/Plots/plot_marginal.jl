const ep_color = :red
const fba_color = :blue
const hr_color = :black

# Single out
function plot_marginal!(p, metnet::MetNet, out::Union{FBAout, EPout}, ider; h = 1,
    color = :black, label = rxns(metnet, ider), kwargs...)

    lb_ = lb(metnet, ider)
    ub_ = ub(metnet, ider)

    dist = marginal(metnet, out, ider)

    if sqrt(var(dist)) == 0.0
        plot!(p, [lb_ - ub_/10, ub_ + ub_/10], [0.0, 0.0]; 
            label = "", color = color, kwargs...)
        plot!(p, [mean(dist), mean(dist)], [0.0, h]; 
            label = label, color = color, kwargs...)
    else
        μ_ = μ(metnet, out, ider);
        σ_ = σ(metnet, out, ider);
        xs = normal_xs_for_plotting(μ_, σ_, lb_ - ub_/10, ub_ + ub_/10)
        ys = [pdf(dist, x) for x in xs]
        plot!(p, xs, ys; label = label, color = color, kwargs...)
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
    plot_marginal!(p, metnet, hrout, ider; color = hr_color, label = "HR")
    plot_marginal!(p, metnet, epout, ider; color = ep_color, lw = 10, label = "EP")
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
    plot_marginal!(p, metnet, epout, ider; color = ep_color, lw = 10, label = "EP")
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

    plot_marginal!(p, metnet, hrout, ider; color = hr_color, label = "HR")
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