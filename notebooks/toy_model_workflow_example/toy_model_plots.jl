# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

using StatsBase
using Serialization
using Plots

# run add "https://github.com/josePereiro/Chemostat" in the julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat;

# ### Loading cache
# This notebook use the cached data from [toy_model_computation.ipynb](./toy_model_computation.ipynb)

cache_file = "toy_model_cache.jls"
bundle = deserialize(cache_file)
sort!(bundle.ξs)
sort!(bundle.βs)
println(cache_file, " loaded!!!")

println("xis: ", bundle.ξs)
println("betas: ", bundle.βs)

# ---
# ## Plots
# --- 

# ### Marginals

# +
# Plot the marginals as function of βs

# An increase β implies that the solutions with greater value of the objective reaction (biom in this case)
# will be more probables (more prob density around the large values of biom). The rest of the fluxes will change 
# depending in its covarianze with biom

ξ = bundle.ξs[4] # Select ξ to plot
metnet = Ch.Utils.get_data(bundle, ξ, :net)
fbaout = Ch.Utils.get_data(bundle, ξ, :fba)
ps = []
βs_ = bundle.βs[6:end] # Select βs to plot
βs_ = [βs_; reverse(βs_)]
iders_ = metnet.rxns

function plot1(iders, ξ, β)
    ps = []
    for ider in iders
        ξstr = round(ξ, digits = 1)
        βstr = round(β, digits = 1)
        p = Plots.plot(title = "$(ider)\nxi: $ξstr, beta: $βstr", 
            xlabel = "flx", ylabel = "pdf", 
            legend = false, yaxis = nothing)
        
        Ch.Plots.plot_marginal!(p, bundle, ξ, β, [:fba, :hr, :ep], ider)
        
        pdf_maxval = Ch.Utils.pdf_maxval(bundle, ξ, β, [:ep, :fba], ider)
        Plots.plot!(p, yaxis = [0.0, pdf_maxval * 1.2])
        push!(ps, p)
    end
    p = Plots.plot(framestyle = :none, 
        legend = :topleft, legendtitlefont = 15, legendtitle = "Toy model", legendfont = 15)
    p = Plots.plot!(p, [],[], lw = 5, color = :blue, label = "FBA")
    Plots.plot!(p, [],[], lw = 5, color = :red, label = "EP")
    Plots.plot!(p, [],[], lw = 5, color = :red, label = "HR")
    
    push!(ps, p);
    
    return Plots.plot(ps..., size = [1000, 1000])
end

gif = Plots.@gif for i in eachindex(βs_)
    plot1(iders_, ξ, βs_[i])
end every 1

# saving gif
gif_file = "toy_model__all_marginals__varying_beta.gif"
cp(gif.filename, gif_file, force = true)
println(relpath(gif_file), " created!!!")
gif
# -

rm(gif.filename, force = true);

# +
# Plot the marginals as function of ξs

# Because increase xi implies a reduction in the intakes all fluxes of the network
# tends to the closest to zero limit possible

β = bundle.βs[5] # Select β to plot
ξs_ = bundle.ξs # Select ξs to plot
ξs_ = [ξs_; reverse(ξs_)]
iders_ = metnet.rxns # Select rxn iders

function plot2(iders, ξ, β, ξ0)
    ps = []
    
    metnet0 = Ch.Utils.get_data(bundle, ξ0, :net)
    
    for ider in iders
        ξstr = round(ξ, digits = 1)
        βstr = round(β, digits = 1)
        
        p = Plots.plot(title = "$(ider)\nxi: $ξstr, beta: $βstr", 
            xlabel = "flx", ylabel = "pdf", 
            legend = false,
            yaxis = nothing)
        
        Ch.Plots.plot_marginal!(p, bundle, ξ, β, [:fba, :hr, :ep], ider)
        
        lb_, ub_ = Ch.Utils.bounds(metnet0, ider)
        margin_ = abs(ub_ - lb_) * 0.1
        pdf_maxval = Ch.Utils.pdf_maxval(bundle, ξ, β, [:ep, :fba], ider)
        Plots.plot!(p, 
            xaxis = [lb_ - margin_, ub_ + margin_],
            yaxis = [0.0, pdf_maxval * 1.2])
        
        push!(ps, p)
    end
    p = Plots.plot(framestyle = :none, 
        legend = :topleft, legendtitlefont = 15, legendtitle = "Toy model", legendfont = 15)
    p = Plots.plot!(p, [],[], lw = 5, color = :blue, label = "FBA")
    Plots.plot!(p, [],[], lw = 5, color = :red, label = "EP")
    Plots.plot!(p, [],[], lw = 5, color = :black, label = "HR")
    
    push!(ps, p);
    
    return Plots.plot(ps..., size = [1000, 1000])
end

# gif
gif = Plots.@gif for i in eachindex(ξs_)
    plot2(iders_, ξs_[i], β, minimum(ξs_))
end every 1

# saving gif
gif_file = "toy_model__all_marginals__varying_xi.gif"
cp(gif.filename, gif_file, force = true)
println(relpath(gif_file), " created!!!")
gif
# -

rm(gif.filename, force = true);

# +
# Plot stoi err asfunctio of xi
β = bundle.βs[1] # Select β to plot

ξs_ = bundle.ξs#[[1, 10, 30, 40, 50]] # Select ξs to plot

p = Plots.plot(title = "Toy Model\nbeta: $β", 
    xlabel = "xi", ylabel = "stoi err/ flx", 
)

iders_ = metnet.mets # Select iders to plot
for ider in iders_
    errs_ = [Ch.Utils.stoi_err(bundle, ξ, β, :ep, ider) for ξ in ξs_]
    Plots.plot!(p, ξs_, errs_, label = "", color= :black)
end
Plots.plot!(p, [], [], color = :black, label = "abs stoi err", lw = 3)

flxs_ = [abs.(Ch.Utils.av(bundle, ξ, β, :ep)) for ξ in ξs_]
Plots.plot!(p, ξs_, mean.(flxs_), lw = 3, label = "", color = :white)
Plots.plot!(p, ξs_, mean.(flxs_), ls = :dash, lw = 3, label = "mean abs flx", color = :blue)
Plots.plot!(p, ξs_, minimum.(flxs_), lw = 3, label = "", color = :white)
Plots.plot!(p, ξs_, minimum.(flxs_), ls = :dash, lw = 3, label = "min abs flx", color = :red)

# gif (making a zoom)
yulims_ = Ch.Utils.logspace(-1,1, 100) |> reverse
yulims_ = [fill(maximum(yulims_), 25); yulims_]
yulims_ = [reverse(yulims_); yulims_]
yllim_ = -0.01
gif = Plots.@gif for i in eachindex(yulims_)
    iter_ = yllim_:((yulims_[i] - yllim_)/8):yulims_[i]
    Plots.plot!(p, yaxis = [yllim_, yulims_[i]], yticks = (iter_, 
            string.(round.(collect(iter_), digits = 2))))
end every 1

# saving gif
gif_file = "toy_model__stoi_err_vs_xi.gif"
cp(gif.filename, gif_file, force = true)
println(relpath(gif_file), " created!!!")
gif
# -
rm(gif.filename, force = true);

# +
# Plot stoi err asfunctio of xi
ξ = bundle.ξs[10] # Select ξ to plot
ξstr = round(ξ, digits = 2)
βs_ = bundle.βs[2:end] # Select βs to plot

p = Plots.plot(title = "Toy Model\nxi: $ξstr", 
    xlabel = "beta", ylabel = "stoi err", 
    xaxis = :log10,
#     yaxis = [-1e-3, 1e-3]
)

iders_ = metnet.mets # Select iders to plot
for ider in iders_
    errs_ = [Ch.Utils.stoi_err(bundle, ξ, β, :ep, ider) for β in βs_]
    Plots.plot!(p, βs_, errs_, color = :black, label = "")
end

Plots.plot!(p, [], [], color = :black, label = "abs stoi err")

flxs_ = [abs.(Ch.Utils.av(bundle, ξ, β, :ep)) for β in βs_]
Plots.plot!(p, βs_, mean.(flxs_), lw = 3, label = "", color = :white)
Plots.plot!(p, βs_, mean.(flxs_), ls = :dash, lw = 3, label = "mean abs flx", color = :blue)
Plots.plot!(p, βs_, minimum.(flxs_), lw = 3, label = "", color = :white)
Plots.plot!(p, βs_, minimum.(flxs_), ls = :dash, lw = 3, label = "min abs flx", color = :red)

# gif
yulims_ = Ch.Utils.logspace(-0.5,0.9, 100) |> reverse
yulims_ = [fill(maximum(yulims_), 25); yulims_]
yulims_ = [reverse(yulims_); yulims_]
yllim_ = -0.01
gif = Plots.@gif for i in eachindex(yulims_)
    iter_ = yllim_:((yulims_[i] - yllim_)/8):yulims_[i]
    Plots.plot!(p, yaxis = [yllim_, yulims_[i]], yticks = (iter_, 
            string.(round.(collect(iter_), digits = 2))))
end every 1

# saving gif
gif_file = "toy_model__stoi_err_vs_beta.gif"
cp(gif.filename, gif_file, force = true)
println(relpath(gif_file), " created!!!")
gif
# -
rm(gif.filename, force = true);


