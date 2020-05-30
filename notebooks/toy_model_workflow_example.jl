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

import Plots
using StatsBase

# run add "https://github.com/josePereiro/Ch" in the julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat;

# ---
# ## Meta
# ---

notebook_name = "toy_model_workflow_example"

# ---
# ## EP, FBA, HR
# ---

obj_ider = "biom";
α = 1e9

intake_info = Dict("gt" => Dict("ub" => 100.0, "c" => 10.0));

ξs = Ch.Utils.logspace(-1,1);
ξs_str = round.(ξs, digits = 2);
println(βs_str)

β0 = 0.0;
βs = [β0; Ch.Utils.logspace(0,3.8, 100)] # This is the working interval of HR
βs_str = round.(βs, digits = 2);
println(βs_str)

# +
# TODO add HR
boundle = Ch.Utils.ChstatBoundle()
# TODO check adding custom info to boundle
#     Ch.Utils.add_data!(boundle, "intake_info", intake_info)

progress = [0, length(ξs) * length(βs)]
@time @sync for (ξi, ξ) in enumerate(ξs)
    @async begin

        # Compute β = 0.0
        # Preparing model
        model = Ch.Utils.toy_model();
        model = Ch.SteadyState.apply_bound!(model, ξ, intake_info)
        model = Ch.Utils.fva_preprocess(model)
        obj_idx = Ch.Utils.rxnindex(model, obj_ider)

        fbaout = Ch.FBA.fba(model, obj_ider)

        Ch.Utils.add_data!(boundle, ξ, :fba, fbaout)
        Ch.Utils.add_data!(boundle, ξ, :net, model)    

        # Compute β > 0.0
        βv = zeros(size(model, 2))
        seed_epout = nothing
        for (βi, β) in enumerate(βs)
            
            βv[obj_idx] = β
            epout = Ch.MaxEntEP.maxent_ep(model, α = α, βv = βv, 
                epsconv = 1e-10, # This is an exageration!!
                solution = seed_epout, verbose = false)
            seed_epout = epout
            
            print("Progress: $(progress[1]/ progress[2])               \r"); flush(stdout);
            progress[1] += 1
            
            Ch.Utils.add_data!(boundle, ξ, β, :ep, epout)
        end
        
    end # @async
end # ξ loop

println()
println("Done!!!", " "^100);
# -

# ---
# ## Plots
# --- 

# ### Marginals

# Plot the marginals as function of β
ξ = boundle.ξs[1] # Select ξ to plot
metnet = Ch.Utils.get_data(boundle, ξ, :net)
fbaout = Ch.Utils.get_data(boundle, ξ, :fba)
ps = []
βs_ = boundle.βs[1:2:end] # Select βs to plot
colors = Plots.distinguishable_colors(length(βs_))
iders_ = metnet.rxns # Select iders to plot
for ider in iders_
    p = Plots.plot(title = ider, 
        legend = false, xlabel = "flx", ylabel = "pdf")

    pdf_maxval = 0.0
    for (i, β) in enumerate(βs_)
        epout = Ch.Utils.get_data(boundle, ξ, β, :ep)

        pdf_maxval = max(pdf_maxval, Ch.Utils.pdf_maxval(metnet, [epout], ider))
        Ch.Plots.plot_marginal!(p, metnet, epout, ider, 
            color = Plots.Gray(1 - i/length(βs_)), 
            alpha = (i/length(βs_)) * 255, 
            label = "", lw = 1)
    end
    Ch.Plots.plot_marginal!(p, metnet, fbaout, ider, color = :blue, label = "FBA", lw = 3, 
        h = pdf_maxval * 1.1, ls = :dash)
    push!(ps, p)
end
# legend
p = Plots.plot(framestyle = :none, legend = :topleft)
p = Plots.plot!(p, [],[], ls = :dash, c = :black, label = "fba")
Plots.plot!(p, [],[], ls = :solid, c = :black, label = "ep")
push!(ps, p);

# An increase β implies that the solutions with greater value of the objective reaction (biom in this case)
# will be more probables (more prob density around the large values of biom). The rest of the fluxes will change 
# depending in its covarianze with biom
Plots.plot(ps..., size = [1000, 1000])

# +
# Plot the marginals as function of ξs
β = boundle.βs[1] # Select β to plot
ps = []
ξs_ = boundle.ξs#[[1, 10, 30, 40, 50]] # Select ξs to plot

iders_ = metnet.rxns # Select iders to plot
for ider in iders_
    p = Plots.plot(title = ider, 
        legend = false, xlabel = "flx", ylabel = "pdf")

    pdf_maxval = 0.0
    for (i, ξ) in enumerate(sort(ξs_))
        metnet = Ch.Utils.get_data(boundle, ξ, :net)
#         fbaout = Ch.Utils.get_data(boundle, ξ, :fba)

        epout = Ch.Utils.get_data(boundle, ξ, β, :ep)
        
        
        pdf_maxval = max(pdf_maxval, Ch.Utils.pdf_maxval(metnet, [epout], ider))
        Ch.Plots.plot_marginal!(p, metnet, epout, ider, 
            color = Plots.Gray(1 - i/length(ξs_)), 
            alpha = (i/length(ξs_)) * 255, 
            label = "", lw = 1)

    end
    push!(ps, p)
end
# legend
p = Plots.plot(framestyle = :none, legend = :topleft)
p = Plots.plot!(p, [],[], ls = :dash, c = :black, label = "fba")
Plots.plot!(p, [],[], ls = :solid, c = :black, label = "ep")
push!(ps, p);
# -

# Because increase xi implies a reduction in the intakes all fluxes of the network
# tends to the closest to zero limit possible
Plots.plot(ps..., size = [1000, 1000])

# +
# Plot stoi err asfunctio of xi
β = boundle.βs[1] # Select β to plot

ξs_ = boundle.ξs#[[1, 10, 30, 40, 50]] # Select ξs to plot

p = Plots.plot(title = "beta: $β", 
    xlabel = "xi", ylabel = "stoi err/ flx", 
#     yaxis = [-1e-4, 1e-4]
)

iders_ = metnet.mets # Select iders to plot
for ider in iders_
    errs_ = [Ch.Utils.stoi_err(boundle, ξ, β, :ep, ider) for ξ in ξs_]
    Plots.plot!(p, ξs_, errs_, label = "", color= :black)
end
Plots.plot!(p, [], [], color = :black, label = "err")

flxs_ = [abs.(Ch.Utils.av(boundle, ξ, β, :ep)) for ξ in ξs_]
Plots.plot!(p, ξs_, minimum.(flxs_), ls = :dash, label = "", color = :black)
Plots.plot!(p, ξs_, mean.(flxs_), ls = :dash, label = "abs av min/mean/max", color = :black)
Plots.plot!(p, ξs_, maximum.(flxs_), ls = :dash, label = "", color = :black)
p

# +
# Plot stoi err asfunctio of xi
ξ = boundle.ξs[10] # Select ξ to plot
βs_ = boundle.βs # Select βs to plotξs_ = boundle.ξs#[[1, 10, 30, 40, 50]] # Select ξs to plot

p = Plots.plot(title = "xi: $ξ", 
    xlabel = "beta", ylabel = "stoi err", 
#     yaxis = [-1e-3, 1e-3]
)

iders_ = metnet.mets # Select iders to plot
for ider in iders_
    errs_ = [Ch.Utils.stoi_err(boundle, ξ, β, :ep, ider) for β in βs_]
    Plots.plot!(p, βs_, errs_, color = :black, label = "")
end

Plots.plot!(p, [], [], color = :black, label = "err")

flxs_ = [abs.(Ch.Utils.av(boundle, ξ, β, :ep)) for β in βs_]
Plots.plot!(p, βs_, minimum.(flxs_), ls = :dash, label = "", color = :black)
Plots.plot!(p, βs_, mean.(flxs_), ls = :dash, label = "abs av min/mean/max", color = :black)
Plots.plot!(p, βs_, maximum.(flxs_), ls = :dash, label = "", color = :black)

p
