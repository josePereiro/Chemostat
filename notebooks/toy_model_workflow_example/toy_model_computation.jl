# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl:light,ipynb
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

import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")
Base.active_project() |> println


# +
using Distributed
using Serialization
using Dates

nworkers = 1 # length(Sys.cpu_info()) - 1
length(workers()) < nworkers && addprocs(nworkers; exeflags = "--project")
println("Working in: ", workers())
# -
@everywhere begin
    using Serialization
    using Distributed
    
    import Chemostat
    Ch = Chemostat
end

# ### Initializing everywhere

@everywhere begin
    
    const obj_ider = "biom";
    const cost_ider ="tot_cost"

    const ep_alpha = 1e9
    const ep_epsconv = 1e-9 # This is an exageration!!
    const ep_ochlen = 100
    
    const hr_b0_nsamples = 100_000 #10_000_000
    const hr_nsamples = 1000 # 10_000
    const intake_info = Dict("gt" => Dict("ub" => 100.0, "c" => 10.0));

    const ξs = Ch.Utils.logspace(-1,1,20);
    const ξs_str = round.(ξs, digits = 2);
    
    const β0 = 0.0 # It is important to include beta = 0 if doing HR
    const βs = [β0; Ch.Utils.logspace(-1, 3.8, 20)]
    const βs_str = round.(βs, digits = 2);
end

println("ep_alpha: ", ep_alpha)
println("ep_epsconv: ", ep_epsconv)
println("βs: ", βs_str)
println("ξs: ", ξs_str)

# RES IDS
# Collect all the computed results ids for bundling
const chnl = RemoteChannel() do
    Channel{Any}(10)
end
const res_ids = []
const collector = @async while true
    id = take!(chnl)
    push!(res_ids, id)
end;

@everywhere begin
    const base_model = Ch.Test.toy_model();
    function prepare_model(ξ)
        model = base_model |> deepcopy
        model = Ch.SteadyState.apply_bound!(model, ξ, intake_info |> deepcopy)
        model = Ch.LP.fva_preprocess(model; verbose = false)
    end
end

# cache dir
cache_dir = "cache"
!isdir(cache_dir) && mkdir(cache_dir)
@assert isdir(cache_dir)
Ch.Utils.set_cache_dir(cache_dir)

# ### EP

# SIMULATION
# Any of the loops can be parallelized by just 
# changing one of the 'map' functions
map(ξs) do ξ

    # HASH SIM
    sim_hash = (:TOY_SIMULATION, ξ)
    
    # FBA EP
    dat = Ch.SimulationUtils.cached_simulation(;
        epochlen = ep_ochlen, 
        verbose = true,
        sim_id = sim_hash,
        get_model = () -> prepare_model(ξ),
        objider = obj_ider, 
        beta_info = [(obj_ider, βs)],
        costider = cost_ider,
        clear_cache = false,
        use_seed = true,
        epmodel_kwargs = Dict(:alpha => ep_alpha),
        epconv_kwargs = Dict(:epsconv => ep_epsconv)
    )
    
    ## HR
    model = prepare_model(ξ)
    hrout_β0 = nothing # seed
    Ch.Utils.tagprintln_inmw("DOING HR ")
    for (βi, β) in enumerate(βs)
        
        Ch.Utils.tagprintln_inmw("DOING BETA ", βi)
        if β == β0 # hrout seed
                hrout_β0 = Ch.MaxEntHR.maxent_hr(model, 
                            nsamples = hr_b0_nsamples, 
                            drop_samples = false, verbose = false);
            hrout = Ch.Utils.HRout(hrout_β0, drop_samples = true) # a copy without samples
        else        
            isnothing(hrout_β0) && (error("The minimum β must be cero!!!"); return)
            
            hrout = Ch.MaxEntHR.maxent_hr(model, hrout_β0, obj_ider, 
                β * 10.0^(-sqrt(ℯ)); # I have no idea why these work! TODO find it!
                nsamples = hr_nsamples, 
                maxiter = 1e16, 
                verbose = false
            )
        end
        dat[(:hr, βi)] = hrout
    end

    ## SAVING DATA
    res_id = (:RESULT, sim_hash)
    Ch.Utils.save_cache(res_id, (ξ, βs, model, dat); 
        headline = "CATCHING RESULTS\n")

    ## PASSING ID TO MASTER
    put!(chnl, res_id)

    GC.gc()
    return nothing
end; # map(ξs) do ξ

for ξ in ξs
    m = prepare_model(ξ)
    Ch.LP.fba(m, obj_ider, cost_ider)
end

# ### Parallel loop

# COLLECTING RESULTS
Ch.Utils.tagprintln_inmw("COLLECTING RESULTS ")
sleep(1) # wait for collector to get all ids
const bundle = Ch.Utils.ChstatBundle()
for id in res_ids

    ξ_, βs_, model_, dat_ = Ch.Utils.load_cache(id; verbose = false)

    bundle[ξ_, :net] = model_
    bundle[ξ_, :fba] = dat_[:fba]

    for (βi, β_) in βs_ |> enumerate
        bundle[ξ_, β_, :ep] = dat_[(:ep, βi)]
        bundle[ξ_, β_, :hr] = dat_[(:hr, βi)]
    end
end

# ### Saving results

### Catching
cache_file = "toy_model_cache.jls"
serialize(cache_file, bundle)
println(relpath(cache_file), " created!!!")

# ### Deleting temporal caches

Ch.Utils.tagprintln_inmw("DELETING TEMPORAL CACHE ")
Ch.Utils.delete_temp_caches()


