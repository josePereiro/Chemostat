# given a model and a new set of bounds check
# how they will affect a given obj reaction
function check_newbounds(S, b, lb, ub, newlb, newub, 
        check_obj::Int, idxs = eachindex(lb); 
        check_obj_atol = 1e-4,
        verbose = true,
        batchlen = 50)

    M, N = size(S)
    T = eltype(S)
    nths = nthreads()
    ref_obj_val = fba(S, b, lb, ub, check_obj).obj_val

    # thread environment (avoid race)
    env_pool = Dict()
    for tid in 1:nths
        env = get!(env_pool, tid, Dict())
        env[:newlb] = Dict(i => val for (i, val) in enumerate(newlb))
        env[:newub] = Dict(i => val for (i, val) in enumerate(newub))
        env[:wlb] = deepcopy(lb)
        env[:wub] = deepcopy(ub)
    end

    icount = length(idxs)
    batchlen = max(1, min(batchlen, length(idxs)))
    batches = [idxs[i0:(min(i0 + batchlen - 1, icount))] for i0 in 1:batchlen:icount]
    verbose && (prog = Progress(icount; desc = "Checking bounds  "))
    @threads for batch in batches
        thid = threadid()
        env = env_pool[thid]

        wlb, wub = env[:wlb], env[:wub]
        tnewlb, tnewub = env[:newlb], env[:newub]

        # Test whole batch
        wlb_bkup, wub_bkup = wlb[batch], wub[batch] # backup
        wlb[batch], wub[batch] = [tnewlb[idx] for idx in batch], [tnewub[idx] for idx in batch]
        new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
        wlb[batch], wub[batch] = wlb_bkup, wub_bkup # reset
        if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)

            # If whole batch fail, test independent idxs
            for idx in batch
                
                @assert all(wlb .== lb) # Dev
                @assert all(wub .== ub) # Dev

                # check both first (this use the premise that only a few rxns will affect the biomass)
                wlb_bkup, wub_bkup = wlb[idx], wub[idx] # backup
                wlb[idx], wub[idx] = tnewlb[idx], tnewub[idx]
                new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
                wlb[idx], wub[idx] = wlb_bkup, wub_bkup # reset

                if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                    # if obj_val changed I check each one 
                    for (wcol, fvacol, bkup) in [(wlb, tnewlb, wlb_bkup), 
                                                (wub, tnewub, wub_bkup)]
                        wcol[idx] = fvacol[idx]
                        new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
                        if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                            fvacol[idx] = bkup # reset
                        end
                        wcol[idx] = bkup # reset  
                    end
                end
                verbose && next!(prog)
            end
        else
            for idx in batch
                verbose && next!(prog)
            end
        end
    end # for batch in batches
    verbose && finish!(prog)

    # Collect results
    merged_fvalb, merged_fvaub = Dict{Int, T}(), Dict{Int, T}()
    for tid in 1:nths
        env = get(env_pool, tid, Dict())
        merge!(merged_fvalb, get(env, :newlb, Dict()))
        merge!(merged_fvaub, get(env, :newub, Dict()))
    end

    return [merged_fvalb[idx] for idx in idxs], [merged_fvaub[idx] for idx in idxs]
end
