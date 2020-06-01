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

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES && addprocs(NO_CORES)

@everywhere begin
    using ParallelDataTransfer
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
    α = 1e9
    intake_info = Dict("gt" => Dict("ub" => 100.0, "c" => 10.0));

    ξs = Ch.Utils.logspace(-1,1,20);
    ξs_str = round.(ξs, digits = 2);
    
    β0 = 0.0 # It is important to include beta = 0
    βs = [β0; Ch.Utils.logspace(-1, 3.8, 20)]
    βs_str = round.(βs, digits = 2);
end

println("α: ", α)
println("βs: ", βs_str)
println("ξs: ", ξs_str)

# ### work function

@everywhere function process_xi(ξ)
    
    println("Doing xi: $ξ"); flush(stdout)
    
    data = Dict()
    
    # Preparing model
    model = Ch.Utils.toy_model();
    model = Ch.SteadyState.apply_bound!(model, ξ, intake_info)
    model = Ch.Utils.fva_preprocess(model)
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)

    # fba
    fbaout = Ch.FBA.fba(model, obj_ider)
    
    # hr seed
    hrout_β0 = nothing
    
    # storing
    data[(ξ, :net)] = model
    data[(ξ, :fba)] = fbaout

    βv = zeros(size(model, 2))
    seed_epout = nothing
    for (βi, β) in enumerate(sort(βs))
        
        if β == β0 # hrout seed
            hrout_β0 = Ch.MaxEntHR.maxent_hr(model, nsamples = 10_000_000, 
                                    drop_samples = false, verbose = false);
            hrout = Ch.Utils.HRout(hrout_β0, drop_samples = true) # a copy without samples
        else        
            isnothing(hrout_β0) && (error("The minimum β must be cero!!!"); return)
            hrout = Ch.MaxEntHR.maxent_hr(model, hrout_β0, obj_ider, 
                β * 10.0^(-sqrt(ℯ)); # I have no idea why these work! TODO find it!
                nsamples = 10_000, maxiter = 1e16, 
                verbose = false)
        end

        # epout
        βv[obj_idx] = β
        epout = Ch.MaxEntEP.maxent_ep(model, α = α, βv = βv, 
            epsconv = 1e-10, # This is an exageration!!
            solution = seed_epout, verbose = false)
        seed_epout = epout
        
        # storing
        data[(ξ, β, :ep)] = epout
        data[(ξ, β, :hr)] = hrout

    end # β loop
    
    return data
end

function boundle_data!(boundle, ξ, βs, data)

    Ch.Utils.add_data!(boundle, ξ, :net, data[(ξ, :net)])
    Ch.Utils.add_data!(boundle, ξ, :fba, data[(ξ, :fba)])
    
    for β in βs
        Ch.Utils.add_data!(boundle, ξ, β, :ep, data[(ξ, β, :ep)])
        Ch.Utils.add_data!(boundle, ξ, β, :hr, data[(ξ, β, :hr)])
    end
end

# ### Parallel loop

# +
# this can take a while!!!
boundle = Ch.Utils.ChstatBoundle()

@sync for ξ in ξs
    @async begin
        data = fetch(@spawn process_xi(ξ))
        boundle_data!(boundle, ξ, βs, data)
        flush(stdout)
    end
end

println("Done!!!", " "^100);
# -

### Catching
cache_file = "toy_model_cache.jls"
serialize(cache_file, boundle)
println(relpath(cache_file), " created!!!")


