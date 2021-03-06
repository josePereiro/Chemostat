# I break down ep in several epoch to be able to cache partial results
function epoch_converge_ep!(epmodel::EPModel; 
        epochlen::Int = 10,
        maxiter::Int=2000,  
        onerr = (epout, err) -> rethrow(err),
        before_epoch = (epout) -> (false, nothing),
        after_epoch = (epout) -> (false, nothing),
        epkwargs...
    )
    
    epout = nothing
    curr_iter = 0

    try 
        while (isnothing(epout) || epout.iter == 0) || # First time
                (maxiter > curr_iter &&                # Till maxiter 
                epout.status != :converged)            # Or converged
            
            # For fast feed back the first epoch has length 1
            epochlen_ = (isnothing(epout) || epout.iter == 0) ? 1 : epochlen

            ret, dat = before_epoch(epout)
            ret && return dat

            epout = converge_ep!(epmodel; epkwargs...,
                        verbose = false,
                        iter0 = curr_iter,
                        maxiter = curr_iter + epochlen_, 
                        drop_epfields = false # We need the epfields to update the epmodel
                    )

            curr_iter = epout.iter
            
            ret, dat = after_epoch(epout)
            ret && return dat

        end # while

    catch err

        err isa InterruptException && rethrow(err)

        ret, dat = onerr(epout, err)
        ret && return dat

    end # try

    return epout
    
end