function fva(S, b, lb, ub, idxs = eachindex(lb); 
        check_obj = nothing, check_obj_atol = 1e-4,
        verbose = true, 
        zeroth = 1e-10,
        on_empty_sol = (idx, sense) -> error("FBA failed, empty solution returned!!!"))

    M, N = size(S)
    T = eltype(S)
    nths = nthreads()

    # pools (avoid race)
    env_pool = Dict()
    for tid in 1:nths
        env = get!(env_pool, tid, Dict())
        env[:sv] = zeros(T, N)
        env[:fvalb] = Dict{Int, T}()
        env[:fvaub] = Dict{Int, T}()
    end
    
    if !isnothing(check_obj)
        ref_obj_val = fba(S, b, lb, ub, check_obj).obj_val
        # working bounds (avoid race)
        for tid in 1:nths
            env = get!(env_pool, tid, Dict())
            env[:wlb] = deepcopy(lb)
            env[:wub] = deepcopy(ub)

            @assert all(env[:wlb] .== lb)
            @assert all(env[:wub] .== ub)
        end
    end
    
    n = length(idxs)
    verbose && (prog = Progress(n; desc = "Doing FVA  "))

    to_map = collect(idxs)
    @threads for idx in to_map
        tid = threadid()
        
        # get thread environment
        env = env_pool[tid]
        fvalb, fvaub, sv = env[:fvalb], env[:fvaub], env[:sv]
         
        verbose && next!(prog)
        for (fvacol, sense) in [(fvalb, 1.0), (fvaub, -1.0)]
            
            sv[idx] = sense
            sol = linprog(
                sv, # Opt sense vector 
                S, # Stoichiometric matrix
                b, # row lb
                b, # row ub
                lb, # column lb
                ub, # column ub
                ClpSolver());
            x = isempty(sol.sol) ? on_empty_sol(idx, sense) : sol.sol[idx]
            fvacol[idx] = abs(x) < zeroth ? zero(x) : x
            sv[idx] = zero(sense)
        end

        if !isnothing(check_obj)

            wlb, wub = env[:wlb], env[:wub]
            @assert all(wlb .== lb)
            @assert all(wub .== ub)

            # check both first (this use the premise that only a few rxns will affect the biomass)
            wlb_bkup, wub_bkup = wlb[idx], wub[idx] # backup
            wlb[idx], wub[idx] = fvalb[idx], fvaub[idx]
            new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
            wlb[idx], wub[idx] = wlb_bkup, wub_bkup # reset
            if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                # if obj_val changed I check each one 
                for (wcol, fvacol, bkup) in [(wlb, fvalb, wlb_bkup), 
                                             (wub, fvaub, wub_bkup)]
                    wcol[idx] = fvacol[idx]
                    new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
                    if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                        fvacol[idx] = bkup # reset
                    end
                    wcol[idx] = bkup # reset  
                end
            end
        end
    end # for idx in idxs

    verbose && finish!(prog)

    # Collect results
    merged_fvalb, merged_fvaub = Dict{Int, T}(), Dict{Int, T}()
    for tid in 1:nths
        env = get(env_pool, tid, Dict())
        merge!(merged_fvalb, get(env, :fvalb, Dict()))
        merge!(merged_fvaub, get(env, :fvaub, Dict()))
    end

    # return just the indexed bounds
    return [merged_fvalb[idx] for idx in idxs], [merged_fvaub[idx] for idx in idxs]
end

function fva(model::MetNet, iders = eachindex(model.lb); 
        check_obj = nothing, kwargs...) 
    obj_idx = isnothing(check_obj) ? nothing : rxnindex(model, check_obj)
    idxs = [rxnindex(model, idx) for idx in iders]
    model_fields = _extract_dense(model, [:S, :b, :lb, :ub])
    return fva(model_fields..., idxs; check_obj = obj_idx, kwargs...);
end

fva(model::MetNet, ider::IDER_TYPE; kwargs...) = fva(model, [ider]; kwargs...)