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

# +
using Distributed
using Serialization
using Dates

# Precompiling in worker 1
using Chemostat
Ch = Chemostat

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES && addprocs(NO_CORES - 1)

@everywhere begin
    using Serialization
    using Chemostat
    Ch = Chemostat
    using Distributed
end
# -

# ---
# ## EP, FBA, HR
# ---

# ### Initializing everywhere

@everywhere begin
    obj_ider = "biom";
    cost_ider ="tot_cost"
    
    ep_alpha = 1e9
    ep_epsconv = 1e-9 # This is an exageration!!
    
    hr_b0_nsamples = 10_000_000
    hr_nsamples = 10_000
    intake_info = Dict("gt" => Dict("ub" => 100.0, "c" => 10.0));

    ξs = Ch.Utils.logspace(-1,1,20);
    ξs_str = round.(ξs, digits = 2);
    
    β0 = 0.0 # It is important to include beta = 0 if doing HR
    βs = [β0; Ch.Utils.logspace(-1, 3.8, 20)]
    βs_str = round.(βs, digits = 2);
end

println("ep_alpha: ", ep_alpha)
println("ep_epsconv: ", ep_epsconv)
println("βs: ", βs_str)
println("ξs: ", ξs_str)

# ### Print functions

@everywhere function print_hello(wid, xi)
    Core.println("Worker $wid starting xi $xi at $(Time(now())) ----------------------------")
    Core.println("\txis:   $ξs")
    Core.println("\tbetas: $βs")
    Core.println()
    flush(stdout);
end

@everywhere function print_progress(wid, ξi, βi, elapsed)
    Core.println("Worker: $wid xi $ξi at $(Time(now())) ----------------------------")
    Core.println("\txi: [$ξi/ $(length(ξs))] beta: [$βi/ $(length(βs))]")
    Core.println("\t  ----------------- --------------")
    Core.println("\txi:                 $(ξs[ξi])")
    Core.println("\tbeta:               $(βs[βi])")
    Core.println("\tep_alpha:           $ep_alpha")
    Core.println("\tep_epsconv:         $ep_epsconv")
    Core.println("\thr_b0_nsamples:     $hr_b0_nsamples")
    Core.println("\thr_nsamples:        $hr_nsamples")
    
    Core.println("\t  ----------------- --------------")
    Core.println("\telapsed time(s):    $elapsed")
    Core.println()
    flush(stdout);
end

@everywhere function print_good_bye(wid, xi, tcache_file)
    Core.println("Worker $wid finished xi $xi at $(Time(now())) ----------------------------")
    Core.println("\tTemp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

@everywhere function print_return_cached(wid, xi, tcache_file)
    Core.println("Worker $wid returns cached xi $xi at $(Time(now())) ----------------------------")
    Core.println("\tTemp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

# ### temp cache file

@everywhere temp_cache_file(xi) = "toy_model_temp_cache_xi_$xi.jls"

# ### work function

@everywhere function process_xi(ξi; upfrec = 5)
    
    ξ = ξs[ξi]
    
    # I will cache temporally the results of this function 
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    tcache_file = temp_cache_file(ξi)

    # If cached return 
    if isfile(tcache_file)
        (ξ, data) =  deserialize(tcache_file)
        # Info
        # Print in worker 1
        remotecall_wait(print_return_cached, 1, myid(), ξi, tcache_file)
        return (ξ, data)
    end
    
    # Info
    # Print in worker 1
    remotecall_wait(print_hello, 1, myid(), ξi)
    
    # Strating time
    t0 = time()
    
    data = Dict()
    
    # Preparing model
    model = Ch.Utils.toy_model();
    model = Ch.SteadyState.apply_bound!(model, ξ, intake_info)
    model = Ch.Utils.fva_preprocess(model)
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)

    # fba
    fbaout = Ch.FBA.fba(model, obj_ider, cost_ider)
    
    # hr seed
    hrout_β0 = nothing
    
    # storing
    data[(ξ, :net)] = model
    data[(ξ, :fba)] = fbaout

    βv = zeros(size(model, 2))
    seed_epout = nothing
    for (βi, β) in enumerate(sort(βs))
        
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
                verbose = false)
        end

        # epout
        βv[obj_idx] = β
        epout = Ch.MaxEntEP.maxent_ep(model, 
            alpha = ep_alpha, 
            beta_vec = βv, 
            epsconv = ep_epsconv, 
            solution = seed_epout, verbose = false)
        seed_epout = epout
        
        # storing
        data[(ξ, β, :ep)] = epout
        data[(ξ, β, :hr)] = hrout
        
        # Info
        # Print in worker 1
        show_progress = βi == 1 || βi == length(βs) || βi % upfrec == 0
        show_progress && remotecall_wait(print_progress, 1, myid(), ξi, βi, time() - t0)

    end # β loop
    
    # Catching
    serialize(tcache_file, (ξ, data))
    
    # Info
    # Print hello in worker 1
    remotecall_wait(print_good_bye, 1, myid(), ξi, tcache_file)
    
    return (ξ, data)
end

function bundle_data!(bundle, ξ, βs, data)

    Ch.Utils.add_data!(bundle, ξ, :net, data[(ξ, :net)])
    Ch.Utils.add_data!(bundle, ξ, :fba, data[(ξ, :fba)])
    
    for β in βs
        Ch.Utils.add_data!(bundle, ξ, β, :ep, data[(ξ, β, :ep)])
        Ch.Utils.add_data!(bundle, ξ, β, :hr, data[(ξ, β, :hr)])
    end
end

# ### Parallel loop

# this can take a while!!!
remote_results = pmap(process_xi, eachindex(ξs));

# ### Saving results

# +
bundle = Ch.Utils.ChstatBundle()

for (ξ, data) in remote_results
    bundle_data!(bundle, ξ, βs, data)
end

println("Done!!!");
# -

### Catching
cache_file = "toy_model_cache.jls"
serialize(cache_file, bundle)
println(relpath(cache_file), " created!!!")

# ### Deleting temporal caches

println("Deleting temporal cache files")
for ξi in eachindex(ξs)
    tcache_file = temp_cache_file(ξi)
    !isfile(tcache_file) && continue
    rm(tcache_file, force = true)
    println(relpath(tcache_file), " deleted!!!")
end



