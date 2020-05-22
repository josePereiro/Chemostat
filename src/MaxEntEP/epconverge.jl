# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function epconverge!(epfield::EPFields, epmat::M, epalg::EPAlg, eponesweep!::T) where {T<:Function,M<:AbstractEPMat}


    @extract epalg : maxiter verbose epsconv

    returnstatus = :unconverged
    iter = 0
    verbose && print(stderr, "Converging with β=$(epalg.beta) maxth=$(epsconv) maxiter=$(maxiter):\n")
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    while iter < maxiter
        iter += 1
        # eponesweep! will be eponesweepT0! or eponesweep depending on beta
        (errav, errvar, errμ, errs) = eponesweep!(epfield, epalg, epmat)
        if (verbose)
            @printf(stderr, 
                "it = %d ɛav = %.2e ɛvar = %.2e ɛμ = %.2e ɛs = %.2e                                   \r", 
                iter, errav, errvar, errμ, errs)
            flush(stderr)
        end

        if max(errav, errvar) < epsconv
            returnstatus = :converged
            break
        end
    end

    if verbose
        print(stderr, "\n")
        flush(stderr)
    end

    return returnstatus
end