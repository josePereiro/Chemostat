const UNCONVERGED_STATUS = :unconverged
const CONVERGED_STATUS = :converged

function converge_ep!(epmodel::EPModel{T};
        verbose::Bool=true,  # output verbosity
        damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
        epsconv::Real=1e-6,  # convergence criterion
        maxiter::Int=2000,   # maximum iteration count
        maxvar::Real=1e50,   # maximum numerical variance
        minvar::Real=1e-50,  # minimum numerical variance
        iter0 = 0,           # the started iteration count
        drop_epfields = false  # if the final EPout object will export the epfields
        
    ) where {T<:Real}

    # Unpack
    epfields = epmodel.epfields
    epmat = epmodel.epmat
    beta_vec = epmodel.beta_vec
    alpha = epmodel.alpha
    stat = epmodel.stat

    stat[:converge_init_time] = time()
    epalg = EPAlg(alpha, beta_vec, 
        minvar, maxvar, epsconv, damp, maxiter, verbose)

    returnstatus = UNCONVERGED_STATUS
    iter = iter0
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    prog = ProgressThresh{typeof(epsconv)}(epsconv; desc =  "EP  ", dt = 0.5)
    max_beta = findmax(beta_vec)
    for iter in iter0:maxiter

        # eponesweep! will be eponesweepT0! or eponesweep depending on alpha
        stat[:elapsed_eponesweep] = @elapsed begin
            errav, errvar, errμ, errs = epmodel.updatealg!(epfields, epalg, epmat, stat)
        end

        max_err = max(errav, errvar, errμ, errs)
        max_err < epsconv && (returnstatus = CONVERGED_STATUS; break)

        if verbose 
            sweep_time = stat[:elapsed_eponesweep]
            inv_time = stat[:elapsed_eponesweep_inv]
            inv_frac = round(inv_time * 100/ sweep_time; digits = 3)

            update!(prog, max_err; showvalues = [
                (:iter, iter),
                (:maxiter, maxiter),
                (:alpha, alpha),
                (:max_beta, max_beta),
                (:sweep_time, sweep_time),
                (:inv_time, string(inv_time, " [", inv_frac, " %]")),
            ])
        end
    end

    verbose && finish!(prog)

    #= Scale back μ, s, av, va of epfields and lub, llb and Y =#
    scaleepfield!(epmodel.scalefact, epfields)
    μ, σ = epfields.μ, epfields.s
    av, va = epfields.av, epfields.va
    if isinf(alpha)
        sperm = sortperm(epmat.usperm)
        μ, σ, av, va = μ[sperm], σ[sperm], av[sperm], va[sperm]
    end
    
    sol = drop_epfields ? nothing : epfields
    return  EPout(μ, σ, av, va, sol, returnstatus, iter)
end