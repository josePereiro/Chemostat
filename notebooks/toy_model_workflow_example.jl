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

# run add "https://github.com/josePereiro/Chemostat" in the julia Pkg REPL for installing the package
import Chemostat

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

ξs = Chemostat.Utils.logspace(-1,1);
ξs_str = round.(ξs, digits = 2);
println(βs_str)

β0 = 0.0;
βs = [β0; Chemostat.Utils.logspace(0,3.8, 100)] # This is the working interval of HR
βs_str = round.(βs, digits = 2);
println(βs_str)

# +
# TODO add HR
boundle = Chemostat.Utils.ChstatBoundle()
# TODO check adding custom info to boundle
#     Chemostat.Utils.add_data!(boundle, "intake_info", intake_info)

progress = [0, length(ξs) * length(βs)]
@time @sync for (ξi, ξ) in enumerate(ξs)
    @async begin

        # Compute β = 0.0
        # Preparing model
        model = Chemostat.Utils.toy_model();
        model = Chemostat.SteadyState.apply_bound!(model, ξ, intake_info)
        model = Chemostat.Utils.fva_preprocess(model)
        obj_idx = Chemostat.Utils.rxnindex(model, obj_ider)

        fbaout = Chemostat.FBA.fba(model, obj_ider)

        Chemostat.Utils.add_data!(boundle, ξ, :fba, fbaout)
        Chemostat.Utils.add_data!(boundle, ξ, :net, model)    

        # Compute β > 0.0
        βv = zeros(size(model, 2))
        seed_epout = nothing
        for (βi, β) in enumerate(βs)
            
            βv[obj_idx] = β
            epout = Chemostat.MaxEntEP.maxent_ep(model, α = α, βv = βv, 
                solution = seed_epout, verbose = false)
            seed_epout = epout
            
            print("Progress: $(progress[1]/ progress[2])               \r"); flush(stdout);
            progress[1] += 1
            
            Chemostat.Utils.add_data!(boundle, ξ, β, :ep, epout)
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
ξ = boundle.ξs[1]
metnet = Chemostat.Utils.get_data(boundle, ξ, :net)
fbaout = Chemostat.Utils.get_data(boundle, ξ, :fba)
ps = []
βs_ = boundle.βs[1:2:end] # Select βs to plot
colors = Plots.distinguishable_colors(length(βs_))
iders_ = metnet.rxns # Select iders to plot
for ider in iders_
    p = Plots.plot(title = ider, 
        legend = false, xlabel = "flx", ylabel = "pdf")

    pdf_maxval = 0.0
    for (i, β) in enumerate(βs_)
        epout = Chemostat.Utils.get_data(boundle, ξ, β, :ep)

        pdf_maxval = max(pdf_maxval, Chemostat.Utils.pdf_maxval(metnet, [epout], ider))
        Chemostat.Plots.plot_marginal!(p, metnet, epout, ider, 
            color = Plots.Gray(1 - i/length(βs_)), 
            alpha = (i/length(βs_)) * 255, 
            label = "", lw = 1)
    end
    Chemostat.Plots.plot_marginal!(p, metnet, fbaout, ider, color = :blue, label = "FBA", lw = 3, 
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
        metnet = Chemostat.Utils.get_data(boundle, ξ, :net)
#         fbaout = Chemostat.Utils.get_data(boundle, ξ, :fba)

        epout = Chemostat.Utils.get_data(boundle, ξ, β, :ep)
        
        
        pdf_maxval = max(pdf_maxval, Chemostat.Utils.pdf_maxval(metnet, [epout], ider))
        Chemostat.Plots.plot_marginal!(p, metnet, epout, ider, 
            color = Plots.Gray(1 - i/length(ξs_)), 
            alpha = (i/length(ξs_)) * 255, 
            label = "", lw = 1)

    end
#     Chemostat.Plots.plot_marginal!(p, metnet, fbaout, ider, color = :blue, label = "FBA", lw = 3, 
#         h = pdf_maxval * 1.1, ls = :dash)
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


