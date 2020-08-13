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

# +
using Plots
pyplot();

import Chemostat
const Ch = Chemostat;
# -

function ecToyModel(;Inf_ = 1000)
    # Kinetic params
    k13 = rand()
    k14 = rand()
    k24 = rand()
    k5f = rand()
    k5b = rand()

    S = [# v1     v2     v3    v4/arm  v41    v42    v5f    v5b    v6         e1     e2     e3     e4

           +1      0     -1      0      0      0      0      0      0         0      0      0      0 ;  # M1
            0     +1     -1      0      0      0     -1     +1      0         0      0      0      0 ;  # M2
            0      0     +1     -1      0      0      0      0      0         0      0      0      0 ;  # M3
            0      0      0      0     +1     +1     +1     -1     -1         0      0      0      0 ;  # M4
            0      0      0     +1     -1     -1      0      0      0         0      0      0      0 ;  # PM1


            0      0   -1/k13    0   -1/k14    0      0      0      0        +1      0      0      0 ;  # E1
            0      0      0      0      0   -1/k24    0      0      0         0     +1      0      0 ;  # E2
            0      0      0      0      0      0   -1/k5f -1/k5b    0         0      0     +1      0 ;  # E3
            0      0      0      0      0      0   -2/k5f -2/k5b    0         0      0      0     +1 ;  # E4
        ]

    mets = ["M1", "M2", "M3", "M4", "PM1", "E1", "E2", "E3", "E4"]
    rxns = ["v1", "v2", "v3", "v4/arm", "v41", "v42", "v5f", "v5b", "v6",  "e1" ,  "e2" ,  "e3" ,  "e4" ]
    ub =   [Inf_, Inf_, Inf_,   Inf_  , Inf_ , Inf_ , Inf_ , Inf_ , Inf_, rand(), rand(), rand(), rand()]
    lb = zeros(size(S, 2));
    b = zeros(size(S, 1));
    model = Ch.Utils.MetNet(S, b, lb, ub, rxns, mets);
    return Ch.LP.fva_preprocess(model, verbose = false);
end

# <img src="./ecToyModel.png" width="400" height="700">

model = ecToyModel();

println("Balance equations")
foreach(model.mets) do met
    Ch.Utils.balance_str(model, met) |> println
end

println("Reaction equations")
foreach(model.rxns) do rxn
    eq = Ch.Utils.rxn_str(model, rxn) 
    b = Ch.Utils.bounds(model, rxn)
    println(rpad(rxn, 10), " eq: ", eq, "   ", collect(b))
end

# ## FBA

obj_ider = "v6"
fbaout = Ch.LP.fba(model, obj_ider);

let p1 = plot(title = "FBA", xlabel = "flx_id", ylabel = "flux_val"), 
    p2 = plot(title = "MODEL", xlabel = "flx_id", ylabel = "UB", yscale = :log10)
    plot([bar!(p1, model.rxns, fbaout.v, label = "", color = :red), 
          bar!(p2, model.rxns, model.ub, label = "", yscale = :log10, color = :blue)]...,
    size = [900, 400])
end

# ## EP

epout = Ch.MaxEntEP.maxent_ep(model, verbose = false);

# ## MC-HR

hrout = Ch.MaxEntHR.maxent_hr(model; nsamples = 500_000, drop_samples = true);

let ps = map(model.rxns) do rxn
        p = plot(legend = false, title = rxn, xlim = [0.0, 1.0])
        Ch.Plots.plot_marginal!(p, model, [hrout, epout, fbaout], rxn)
    end
    ncols = 4;
    while mod(length(ps), ncols) != 0
        empty_p = plot([],[],title = "", grid = false, label = "", axis = false)
        push!(ps, empty_p)
    end
    n = length(ps)
    nrows = n / ncols |> Int
    plot(ps...; size = [ncols * 400, nrows * 400], 
        layout = grid(nrows, ncols))
end


