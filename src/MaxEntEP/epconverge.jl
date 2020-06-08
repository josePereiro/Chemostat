# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function epconverge!(epfield::EPFields, epmat::M, epalg::EPAlg, eponesweep!::T) where {T<:Function,M<:AbstractEPMat}
    
    @extract epalg : maxiter verbose epsconv

    returnstatus = :unconverged
    iter = 0
    verbose && print(stderr, "Converging with β=$(epalg.alpha) epsconv=$(epsconv) maxiter=$(maxiter):\n")
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    while iter < maxiter
        iter += 1
        # eponesweep! will be eponesweepT0! or eponesweep depending on alpha
        (errav, errvar, errμ, errs) = eponesweep!(epfield, epalg, epmat)
        if (verbose)
            @printf(stderr, 
                "it = %4d ɛav = %4.2e ɛvar = %4.2e      \r", 
                iter, errav, errvar)
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

    return returnstatus, iter
end