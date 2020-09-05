# Perform fba and ep simulations to a given model.
# The method try to have in ram only what is strictly 
# necessary. It will catch intermedite results and detach
# all unnecessary variables. It is important to feed the 
# method with data that is not linked outside of it 
# to really effectively save memory. 
# In the case of EP, partial results will be cached
# with a frequency regulated by 'epochlen'
function cached_simulation(;
        get_model::Function,
        objider,
        beta_info::Vector = [(objider, [0.0])],
        costider = nothing,
        epochlen::Int = 10,
        epmodel_kwargs::Dict = Dict(),
        epconv_kwargs::Dict = Dict(), # kwargs for the ep convergence alg
        fba_kwargs::Dict = Dict(),
        sim_id = 1, # this should uniquely identify the simulation
        on_hello::Function = () -> nothing,
        verbose = true,
        testing = false,
        clear_cache = true,
        cache_dir = nothing,
        use_seed = true
    )

    # ------------------------------------------------------------------
    # CACHE DIR
    !isnothing(cache_dir) && set_cache_dir(cache_dir)

    # ------------------------------------------------------------------
    # TOP LEVEL CACHE
    cached_data = load_cache(sim_id; verbose = verbose)
    !isnothing(cached_data) && return cached_data

    verbose && tagprintln_inmw("STARTING SIMULATION\n", "sim id: ", sim_id, "\n")

    # ------------------------------------------------------------------
    # MODEL  
    model = get_model()

    T = eltype(model.S)
    M, N = size(model)
    
    # ------------------------------------------------------------------
    # FBA
    fbaout = isnothing(costider) ? fba(model, objider) : fba(model, objider, costider)
    objidx = rxnindex(model, objider)
    fba_objval = av(model, fbaout, objider) 
    verbose && tagprintln_inmw("FBA FINISHED", 
        "\nsim id:          ", sim_id, 
        "\nobjider:         ", objider, 
        "\nfba objval:      ", fba_objval,
        "\ncostider:        ", isnothing(costider) ? "Not set" : costider,
        "\nfba costval:     ", isnothing(costider) ? "Not set" : av(model, fbaout, costider),
        "\n"
    )

    # ------------------------------------------------------------------
    # SAVE AND DROP FBAOUT
    fbaout_id = (:FBAOUT_ID, sim_id)
    save_cache(fbaout_id, fbaout; 
        verbose = verbose, headline = "FBAOUT SAVED")
    fbaout = nothing
    GC.gc()

    # ------------------------------------------------------------------
    # EPMODEL
    # This is the main memory consumer structure. It will be use as a
    # container for all the EP computation. 
    epmodel = EPModel(model; epmodel_kwargs...)

    # ------------------------------------------------------------------
    # SAVE EPFIELD SEED
    # If use_seed = false, the epmodel will be reseted using the original 
    # epfields
    epfield_seed_id = (:EFIRLD_SEED, sim_id)
    if !use_seed
        save_cache(epfield_seed_id, epmodel.epfields; verbose = verbose, 
            headline = "EPFIELD SEED SAVED")
    end

    # ------------------------------------------------------------------
    # PROCESS BETA_INFO
    bcount = _beta_count(beta_info)
    beta_info = _indexed_beta_info(model, beta_info) # Change iders for indexes

    # ------------------------------------------------------------------
    # SAVE AND DROP MODEL
    # EP could be memory consuming, so we drop what is not required
    model_id = (:MODEL, sim_id)
    save_cache(model_id, model |> _compress; 
        verbose = verbose, headline = "MODEL SAVED")
    model = nothing
    GC.gc()

    # SOME INFO
    t0 = time()

    # ------------------------------------------------------------------
    # CONVERGE EP
    beta_vec = spzeros(T, N)
    epout_seed_id = (:EPOUT_SEED, sim_id)
    epout_seed_cfile = temp_cache_file(epout_seed_id)
    
    for βi in 1:bcount

        # ------------------------------------------------------------------
        # UPDATE BETA VEC
        # This just add the desire betas to the vector
        _update_beta_vec!(βi, beta_vec, beta_info)

        # ------------------------------------------------------------------
        # TOP LEVEL BETA CACHE
        # It contains a epout. It will be loaded on demand
        top_beta_id = (:TOP, hash(beta_vec), sim_id)
        cfile = temp_cache_file(top_beta_id)
        isfile(cfile) && continue

        # ------------------------------------------------------------------
        # EPOUT SEED CACHE OR PARTIAL BETA CACHE
        # If use_seed = true a epout seed will be use.
        # This allows to start from a known solution, presumably 
        # closer to the nest than a random one
        # if a partial result is available, it will be used
        # to update the epmodel, if not, I'll try to find a 
        # top seed (The result of the first beta). 
        partial_beta_id = (:PARTIAL, hash(beta_vec), sim_id)
        partial_beta_cfile = temp_cache_file(partial_beta_id)

        if βi != 1 # Not the first time
            seed = nothing
            use_epout_seed = use_seed && (isfile(partial_beta_cfile) || isfile(epout_seed_cfile))
            if use_epout_seed
                epout_id = isfile(partial_beta_cfile) ? partial_beta_id : epout_seed_id
                seed = load_cache(epout_id; verbose = verbose, 
                    headline = "EPOUT SEED LOADED")
            end
            
            # I will reset if explicitly say so
            reset_epmodel = !use_seed 
            if reset_epmodel
                seed = load_cache(epfield_seed_id; verbose = verbose, 
                    headline = "EPFIELD SEED LOADED")
            end

            isnothing(seed) && error("Seed not found!!!")
            update_solution!(epmodel, seed)
            seed = nothing
            GC.gc()
        end

        # ------------------------------------------------------------------
        # PROCESSING BETA
        verbose && tagprintln_inmw("STARTING BETA PROCESSING", 
            "\nsim id:      ", sim_id, 
            "\nbeta count:  [", βi, "/", bcount, "]", 
            "\nmax beta:    ", maximum(beta_vec), 
            "\n")

        cur_epout_cfile = nothing
        epout = epoch_converge_ep!(epmodel;
                epochlen = epochlen,

                before_epoch = function(epout)
                    return (false, nothing)
                end,

                after_epoch = function(epout)

                    # ------------------------------------------------------------------
                    # SAVE PARTIAL RESULT
                    save_cache(partial_beta_id, epout; verbose = verbose,
                        headline = "PARTIAL RESULT SAVED")

                    # Computing stoi error
                    model = load_cache(model_id; verbose = false);
                    isnothing(model) && error(string("Model cache is missing, model_id: ", model_id))
                    stoierr = norm_abs_stoi_err(model[:S], epout.av, model[:b])
                    ep_objval = epout.av[objidx]

                    verbose && tagprintln_inmw("EPOCH FINISHED", 
                        "\nsim id:                   ", sim_id, 
                        "\niter:                     ", epout.iter,
                        "\nmax beta:                 ", maximum(beta_vec), 
                        "\nstoi err (min/mean/max):  ", (minimum(stoierr), mean(stoierr), maximum(stoierr)),
                        "\nfba_objval:               ", fba_objval,
                        "\nep_objval:                ", ep_objval,
                        "\nep status:                ", epout.status,
                        "\n")
                    return (false, nothing)
                end

        ) # epoch_converge_ep!

        # ------------------------------------------------------------------
        # CACHE SEED
        if use_seed && βi == 1
            save_cache(epout_seed_id, epout; verbose = verbose, 
                headline = "EPOUT SEED SAVED")
        end

        # ------------------------------------------------------------------
        # CACHE TOP BETA RESULT
        save_cache(top_beta_id, epout; verbose = verbose, 
                headline = "TOP BETA EPOUT SAVED")

    end

    # ------------------------------------------------------------------
    # DROP EPMODEL
    epmodel = nothing
    GC.gc()

    # ------------------------------------------------------------------
    # LOAD RESULTS
    clear_cache && verbose && tagprintln_inmw("RELOADING RESULTS\n")
    simdata = Dict()
    simdata[:fba] = load_cache(fbaout_id; verbose = false)

    for βi in 1:bcount

        # UPDATE BETA VEC
        # I used the beta_vec hash as id
        _update_beta_vec!(βi, beta_vec, beta_info)

        top_beta_id = (:TOP, hash(beta_vec), sim_id)
        simdata[(:ep, βi)] = load_cache(top_beta_id; verbose = false)
    end

    # ------------------------------------------------------------------
    # CLEAR CACHES
    clear_cache && verbose && tagprintln_inmw("CLEARING CACHES", 
        "\ncache_dir: ", cache_dir
    )
    clear_cache && delete_temp_caches(verbose = verbose)

    return simdata
end

_compress(obj) = obj |> struct_to_dict |> compress_dict


## Handling beta vec
function _beta_count(beta_info)
    beta_vecs = map(last, beta_info)
    l = length(first(beta_vecs))
    !all(map(length, beta_vecs) .== l) && error("Not all beta vec has the same length")
    return l
end

function _indexed_beta_info(model, beta_info)
    new_beta_info = []
    for (ider, betas) in beta_info
        push!(new_beta_info, (rxnindex(model, ider), betas))
    end
    return new_beta_info
end

function _update_beta_vec!(i, beta_vec, beta_info)
    for (idx, betas) in beta_info
        beta_vec[idx] = betas[i]
    end
    return beta_vec
end