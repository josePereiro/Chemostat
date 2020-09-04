function converge_ep!(epmodel::EPModel{T};
        verbose::Bool=true,  # output verbosity
        damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
        epsconv::Real=1e-6,  # convergence criterion
        maxiter::Int=2000,   # maximum iteration count
        maxvar::Real=1e50,   # maximum numerical variance
        minvar::Real=1e-50,  # minimum numerical variance
    ) where {T<:Real}

    @extract epmodel : scalefact updatealg! epfield epmat beta_vec alpha

    epalg = EPAlg(alpha, beta_vec, minvar, maxvar, epsconv, damp, maxiter, verbose)

    #= Here is were all the work is done, this function will 
    call updatealg till convergence or maxiter is reached =#
    # returnstatus, iter = epconverge!(epfield, epmat, epalg, updatealg)
    @extract epalg : maxiter verbose epsconv

    returnstatus = :unconverged
    iter = 0
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    # alpha = epalg.alpha epsconv maxiter
    prog = ProgressThresh{typeof(epsconv)}(epsconv; desc =  "EP  ")
    max_beta = findmax(beta_vec)
    while iter < maxiter
        iter += 1
        # eponesweep! will be eponesweepT0! or eponesweep depending on alpha
        (errav, errvar, errμ, errs) = updatealg!(epfield, epalg, epmat)

        max_err = max(errav, errvar)
        if max_err < epsconv
            returnstatus = :converged
            break
        end

        verbose && update!(prog, max_err; showvalues = [
            (:iter, iter),
            (:maxiter, maxiter),
            (:alpha, alpha),
            (:max_beta, max_beta)
        ])
    end

    verbose && (finish!(prog); flush(stderr))

    #= Scale back μ, s, av, va of epfield and lub, llb and Y =#
    scaleepfield!(epfield, scalefact)
    if alpha < Inf
        return  EPout(epfield.μ, epfield.s, epfield.av, epfield.va, epfield, returnstatus, iter)
    else
        idx = epmat.idx
        return  EPout(epfield.μ[idx],epfield.s[idx],epfield.av[idx],epfield.va[idx], epfield, returnstatus, iter)
    end
end