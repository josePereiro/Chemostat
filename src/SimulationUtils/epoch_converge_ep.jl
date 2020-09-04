
function epoch_converge_ep!(epmodel::EPModel; 
        epochlen::Int = 10,
        maxiter::Int=2000,  
        onerr = (err) -> rethrow(err),
        before_epoch = (epout) -> (false, nothing),
        after_epoch = (epout) -> (false, nothing),
        kwargs...
    )

    # --------------------  MAXENT-EP WHILE LOOP  --------------------  
    # I break down ep in several epoch to be able to cache partial results
    epout = nothing
    curr_iter = 0

    try

        while (isnothing(epout) || epout.iter == 0) || # First time
            (maxiter > curr_iter && # Till maxiter 
                epout.status != :converged) # Or converged

            ret, dat = before_epoch(epout)
            ret && return dat

            epout = converge_ep!(epmodel; kwargs...,
                        verbose = false,
                        iter0 = curr_iter,
                        maxiter = curr_iter + epochlen, 
                    )

            curr_iter = epout.iter
            
            ret, dat = after_epoch(epout)
            ret && return dat

        end # while

    catch err

        err isa InterruptException && rethrow(err)

        ret, dat = onerr(err)
        ret && return dat

    end # try

    return epout
    
end