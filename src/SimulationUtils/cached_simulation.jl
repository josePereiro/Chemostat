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
        verbose = true,
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
    # SIMDATA
    simdata = Dict{Any, Any}()
    
    # ------------------------------------------------------------------
    # FBA
    fbaout = isnothing(costider) ? 
        fba(model, objider; fba_kwargs...) : 
        fba(model, objider, costider; fba_kwargs...)
    objidx = rxnindex(model, objider)
    fba_objval = av(model, fbaout, objider) 
    verbose && tagprintln_inmw("FBA FINISHED", 
        "\nsim id:          ", sim_id, 
        "\nmodel size:      ", (M, N),
        "\nobjider:         ", objider, 
        "\nfba objval:      ", fba_objval,
        "\ncostider:        ", isnothing(costider) ? "Not set" : costider,
        "\nfba costval:     ", isnothing(costider) ? "Not set" : av(model, fbaout, costider),
        "\n"
    )

    # ------------------------------------------------------------------
    # SAVE AND DROP FBAOUT
    fbaout_id = (:FBAOUT_ID, sim_id)
    simdata[:fba] = fbaout_id
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
        epmodel.beta_vec .= beta_vec #TODO: package this

        # ------------------------------------------------------------------
        # TOP LEVEL BETA CACHE
        # It contains a epout. It will be loaded on demand
        top_epout_id = (:TOP_EPOUT, hash(beta_vec), sim_id)
        simdata[(:ep, βi)] = top_epout_id
        top_epout_cfile = temp_cache_file(top_epout_id)
        isfile(top_epout_cfile) && continue

        # ------------------------------------------------------------------
        # EPOUT SEED CACHE OR PARTIAL BETA CACHE
        # If use_seed = true a epout seed will be use.
        # This allows to start from a known solution, presumably 
        # closer to the next than a random one.
        # If a partial result is available, it will be used
        # to update the epmodel, if not, I'll try to find a 
        # top seed (The result of the first beta). 
        partial_epout_id = (:PARTIAL_EPOUT, hash(beta_vec), sim_id)
        partial_epout_cfile = temp_cache_file(partial_epout_id)
        best_epout_id = (:BEST_EPOUT, partial_epout_id)
        best_epout_cfile = temp_cache_file(best_epout_id)
        best_err_id = (:BEST_ERR, partial_epout_id)
        best_err_cfile = temp_cache_file(best_err_id)


        # Select the available epout seed
        seed = nothing
        seed_id = isfile(partial_epout_cfile) ? partial_epout_id : 
                    isfile(best_epout_cfile) ? best_epout_id : 
                        isfile(epout_seed_cfile) ? epout_seed_id : 
                            nothing
        
        use_epout_seed = use_seed && !isnothing(seed_id)
        if use_epout_seed
            seed = load_cache(seed_id; verbose = verbose, 
                headline = "EPOUT SEED LOADED")
        end

        # I will reset if explicitly say so
        reset_epmodel = !use_seed 
        if reset_epmodel
            seed = load_cache(epfield_seed_id; verbose = verbose, 
                headline = "EPFIELD SEED LOADED")
        end

        if !isnothing(seed) 
            update_solution!(epmodel, seed)
            seed = nothing
            GC.gc()
        else
            verbose && tagprintln_inmw("NOT SEEDING", 
                "\nsim id: ", sim_id, "\n")
        end

        # ------------------------------------------------------------------
        # PROCESSING BETA
        verbose && tagprintln_inmw("STARTING BETA PROCESSING", 
            "\nsim id:      ", sim_id, 
            "\nmodel size:  ", (M, N),
            "\nbeta count:  [", βi, "/", bcount, "]", 
            "\nmax beta:    ", maximum(beta_vec), 
            "\n")

        cur_epout_cfile = nothing
        epout = epoch_converge_ep!(epmodel;
                epconv_kwargs...,
                epochlen = epochlen,

                before_epoch = function(epout)
                    return (false, nothing)
                end,

                after_epoch = function(epout)

                    # ------------------------------------------------------------------
                    # SAVE PARTIAL RESULT
                    save_cache(partial_epout_id, epout::EPout; verbose = verbose,
                        headline = "PARTIAL RESULT SAVED")

                    # ------------------------------------------------------------------
                    # COMPUTING STOI ERR
                    model = load_cache(model_id; verbose = false);
                    isnothing(model) && error(string("Model cache is missing, model_id: ", model_id))
                    stoierr = norm1_stoi_err(model[:S], epout.av, model[:b])
                    min_stoierr, mean_stoierr, max_stoierr = minimum(stoierr), mean(stoierr), maximum(stoierr)
                    ep_objval = epout.av[objidx]

                    # ------------------------------------------------------------------
                    # SAVE BEST ERR
                    best_err = (isfile(best_err_cfile) && isfile(best_epout_cfile)) ?
                         load_cache(best_err_id; verbose = false) : Inf
                    if max_stoierr < best_err
                        save_cache(best_epout_id, epout::EPout; verbose = verbose,
                            headline = "BEST RESULT SAVED")
                        save_cache(best_err_id, max_stoierr::Number; verbose = false)
                    end

                    # ------------------------------------------------------------------
                    # stats
                    stat = epmodel.stat
                    sweep_time = stat[:elapsed_eponesweep]
                    inv_time = stat[:elapsed_eponesweep_inv]
                    inv_frac = round(inv_time * 100/ sweep_time; digits = 3)

                    verbose && tagprintln_inmw("EPOCH FINISHED", 
                        "\nsim id:                      ", sim_id, 
                        "\nmodel size:                  ", (M, N),
                        "\nep iter:                     ", epout.iter,
                        "\nmax beta:                    ", maximum(beta_vec), 
                        "\nep stoi err (min/mean/max):  ", (min_stoierr, mean_stoierr, max_stoierr),
                        "\nfba_objval:                  ", fba_objval,
                        "\nep_objval:                   ", ep_objval,
                        "\nep status:                   ", epout.status,
                        "\nlast sweep time(s):          ", sweep_time,
                        "\nlast inv! time(s):            ", inv_time, " [", inv_frac, " % of the sweep time]",
                        "\n")
                    return (false, nothing)
                end,

                onerr = function(epout, err)
                    tagprintln_inmw("ERROR DURING EP\n", err_str(err))
                    !isfile(best_epout_cfile) && error("EP could not recover from error. Best epout file is missing!!")
                    return true, load_cache(best_epout_id, verbose = verbose;
                            headline = "BEST EPOUT LOADED AND RETURNED")
                end

        ) # epoch_converge_ep!

        # ------------------------------------------------------------------
        # SEED RESULT
        if use_seed
            save_cache(epout_seed_id, epout::EPout; verbose = verbose, 
                    headline = "EPOUT SEED SAVED")
        end
        # ------------------------------------------------------------------
        # CACHE TOP BETA RESULT
        save_cache(top_epout_id, epout::EPout; verbose = verbose, 
                headline = "TOP BETA EPOUT SAVED")

        # ------------------------------------------------------------------
        # FINISHING
        # Letting a gap in printing for easy find
        verbose && tagprintln_inmw("BETA PROCESSING FINISHED", 
            "\nsim id:          ", sim_id, 
            "\n\n\n\n\n\n\n\n\n\n"
        )

    end # for βi in 1:bcount

    # ------------------------------------------------------------------
    # DROP EPMODEL
    epmodel = nothing
    GC.gc()

    # ------------------------------------------------------------------
    # LOAD RESULTS
    verbose && tagprintln_inmw("RELOADING RESULTS", 
        "\nCount: ", length(simdata), "\n")
    for (k, dat_id) in simdata
        simdata[k] = load_cache(dat_id; verbose = false)
    end

    # ------------------------------------------------------------------
    # CLEAR CACHES
    clear_cache && verbose && tagprintln_inmw("CLEARING CACHES", 
        "\ncache_dir: ", cache_dir
    )
    clear_cache && delete_temp_caches(verbose = verbose)

    return simdata
end

_compress(obj) = obj |> struct_to_dict |> compressed_copy


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