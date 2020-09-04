# Perform fba and ep simulations to a given model.
# The method try to have in ram only what is strictly 
# necessary. It will catch intermediary results and detach
# all unnecessary variables. It is important to feed the 
# method with data that is not linked outside of it 
# to really effectively save memory. 
# In the case of EP, partial results will be cached
# with a frequency regulated by 'ep_epoch'
function cached_simulation(;
        get_model::Function,
        objider,
        costider = nothing,
        ep_epoch::Int = 10,
        ep_kwargs::Dict = Dict(),
        fba_kwargs::Dict = Dict(),
        sim_id = 1, # this should uniquely identify the simulation
        on_hello::Function = () -> nothing,
        verbose = true,
        testing = false
    )

    # TOP LEVEL CACHE
    cached_data = load_cache(sim_id; verbose = verbose)
    !isnothing(cached_data) && return cached_data

    tagprintln_inmw("STARTING SIMULATION\n", "sim id: ", sim_id)

    # MODEL  
    model = get_model()

    # OBJ_IDERS
    objidx = rxnindex(model, objider)
    costidx = isnothing(costider) ? nothing : rxnindex(model, costider)



end


