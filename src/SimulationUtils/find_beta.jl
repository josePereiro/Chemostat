_set_if_something!(container::Nothing, key, value)::EPout = value
_set_if_something!(container::Dict, key, value)::EPout = container[key] = value
_get_if_something(f, container::Nothing, key)::EPout = f()
_get_if_something(f, container::Dict, key)::EPout = get(container, key, f())

function find_beta(model::MetNet; 
        obj_ider::IDER_TYPE, target_objval::Real,
        beta0::Real, beta1::Real,
        verbose = true,
        epouts::Union{Nothing, Dict} = nothing, # a container/cache for the computed epouts
        maxiters::Real = 100,   
        errorth::Real = 0.01,

        # call backs
        after_maxent_ep::Function = (epout) -> nothing,

        # maxent_ep kwargs
        alpha::Real=1e7, damp::Real=0.9,
        epsconv::Real=1e-6, ep_maxiter::Int=2000,   
        maxvar::Real=1e50, minvar::Real=1e-50,  
        expval=nothing
    )

    M, N = size(model)
    obj_idx = rxnindex(model, obj_ider)
    beta_vec = zeros(N)
    maxent_ep_kwarg = (;verbose, maxiter = ep_maxiter, alpha, damp, epsconv, maxvar, minvar, expval)
    
    # Computing sense
    verbose && println("Computing seed")
    beta_vec[obj_idx] = beta0
    epout0 = _get_if_something(epouts, beta0) do
        epout_ = maxent_ep(model; beta_vec, maxent_ep_kwarg...)
        _set_if_something!(epouts, beta0, epout_)
    end
    after_maxent_ep(epout0)
    beta_vec[obj_idx] = beta1
    epout1 = _get_if_something(epouts, beta1) do
        epout_ = maxent_ep(model; solution = epout0, beta_vec, maxent_ep_kwarg...)
        _set_if_something!(epouts, beta1, epout_)
    end
    after_maxent_ep(epout1)
    objval0 = av(model, epout0, obj_idx)
    objval1 = av(model, epout1, obj_idx)
    sense = objval0 < objval1 ? 1.0 : -1.0
    min_objvar = min(objval0, objval1)
    max_objvar = max(objval0, objval1)
    if !(min_objvar < target_objval < max_objvar)
        error("The interval do not include the target_abjval = $(target_objval), " *
                "current interval goes [$(min_objvar), $(max_objvar)]")
    end

    epout_seed = epout0

    verbose && (prog = ProgressThresh(errorth, "Searching beta  "))
    for it in 1:maxiters
        beta1 == beta0 && error("beta1 == beta0, try doing a manual search or relax errorth ($errorth)!")

        beta = abs(beta0 - beta1)/2 + beta0
        beta_vec[obj_idx] = beta
        epout = _get_if_something(epouts, beta) do
            epout_ = maxent_ep(model; beta_vec, maxent_ep_kwarg..., 
                verbose = false, 
                solution = epout_seed
            )
            _set_if_something!(epouts, beta, epout_)
        end
        after_maxent_ep(epout)
        curr_objval = av(model, epout, obj_idx)
        epout_seed = epout

        # Assume a constant monotony
        if sense > 0.0
            if curr_objval < target_objval; beta0 = beta
                else; beta1 = beta 
            end
        else
            if curr_objval < target_objval; beta1 = beta
                else; beta0 = beta 
            end
        end
        
        curr_diff = abs(target_objval - curr_objval)
        curr_error = curr_diff / target_objval
        if verbose
            update!(prog, curr_error; 
                showvalues = [
                    ("sense           ", sense),
                    ("curr_diff       ", curr_diff),
                    ("beta0           ", beta0),
                    ("beta            ", beta),
                    ("beta1           ", beta1),
                    ("objval0         ", objval0),
                    ("target_objval   ", target_objval),
                    ("objval1         ", objval1),
                    ("curr_objval     ", curr_objval),
                ]
            )
        end
        if curr_error < errorth 
            verbose && finish!(prog)
            return beta
        end
    end
    verbose && finish!(prog)
    return nothing
end