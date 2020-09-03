# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

"""
res=metabolicEP(S,b,lb,ub,...)


The output in res is of type `EPout`: there are several fields:
-   `μ::Vector`: A parameter linked to the mean of the posterior probability
-   `σ::Vector`: A parameter linked to the std  of the posterior probability
-   `av::Vector`: The mean posterior probability
-   `va::Vector`: The variance of the posterior probability
-   `sol::EPFields`: The internal field status. From this value we can restart the sampling from a specific state.
-   `status::Symbol`: either ``:converged`` or ``:unconverged``.

Input (required)
----
- `S`: MxN matrix (either sparse or dense) please note that if you input a dense version, the algorithm is slighlty more efficient. Dense matrices can be create from sparse ones with `Matrix(S)`.
- `b`: a vector of M intakes/uptakes
- `lb`: a vector of lengh N of lower bounds.
- `ub`: a vector of lengh N of upper bounds.

Input (optional arguments).
----
- `alpha` (inverse temperature::``Real``): default 10^7;  the zero temperature algorithm is run setting ``beta=Inf``.
- `beta_vec` (MaxEnt distribution inverse temperature vector:: ``Vector{T}``): default ``T[]``
- `verbose` (``true`` or ``false``): default ``true``
- `damp` (∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield): default 0.9
- `epsconv` (convergence criterion): default 1e-6
- `maxiter` (maximum number of iterations): default 2000
- `maxvar`  (threshold on maximum variance): default 1e50
- `minvar`  (threshold on minimum variance): default 1e-50
- `solution` (start from solution. Is of type ``EPout``): default: ``nothing``
- `expval` (fix to posterior probability of mean and/or variance to values): default ``nothing``. expval can be either at ``Tuple{Float64,Float64,Int}`` or a ``Vector{Tuple{Float64,Float64,Int}}``. Values can be fixed as``expval=(0.2,0.4,4)`` meaning that for flux index 4 the mean is set to 0.2 and the variance to 0.4. Fixing more values ``expval=[(0.2, 0.3, 4), (0.4, nothing, 5)]``: in this case, we fix the posterior of flux 4 to 0.2 (mean) and 0.3 (variance), while for flux 5 we fix the mean to 0.4 and we keep the variance free.
"""

function maxent_ep(S::AbstractArray{T,2}, b::Array{T,1}, lb::Array{T,1}, ub::Array{T,1};
        alpha::Real=1e7,      # inverse temperature
        beta_vec::Vector{T} = T[], # maxent inverse temperature vector
        verbose::Bool=true,  # output verbosity
        damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
        epsconv::Real=1e-6,  # convergence criterion
        maxiter::Int=2000,   # maximum iteration count
        maxvar::Real=1e50,   # maximum numerical variance
        minvar::Real=1e-50,  # minimum numerical variance
        solution::Union{EPout{T}, Nothing}=nothing,  # start from a solution
        expval=nothing # fix posterior probability experimental values for std and mean
    ) where T<:Real


    # making  a local copy to rescale
    llb = copy(lb) 
    lub = copy(ub)

    #= 
    Create EPField, set scale factor and return updatealg (updatefunction),
    depending on alpha (eponesweepT0! or eponesweep).
    The scalefactor is just the maximum absolute bound (lb or ub).
    If a solution is not given, the EPfield will be fresh
    =#
    updatealg, scalefact, epfield = prepareinput(S,b,llb,lub, alpha, verbose, solution, expval)

    #=
    Scale down μ, s, av, va of epfield and lub, llb and Y using the 
    previous computed scalefactor.
    If epfield is fresh, it only will have effect on lub, llb and Y
    =#
    scaleepfield!(epfield,lub,llb,b,1.0/scalefact) # rescaling fields in [0,1]

    epmat = alpha < Inf ? EPMat(S, b, llb, lub, alpha) : EPMatT0(S,b,llb, lub)
    beta_vec = prepare_βv(epmat, beta_vec)
    epalg = EPAlg(alpha, beta_vec, minvar, maxvar, epsconv, damp, maxiter, verbose)

    #= Here is were all the work is done, this function will 
    call updatealg till convergence or maxiter is reached =#
    returnstatus, iter = epconverge!(epfield, epmat, epalg, updatealg)

    #= Scale back μ, s, av, va of epfield and lub, llb and Y =#
    scaleepfield!(epfield,lub,llb,b,scalefact)
    if alpha < Inf
        return  EPout(epfield.μ, epfield.s, epfield.av, epfield.va, epfield, returnstatus, iter)
    else
        idx = epmat.idx
        return  EPout(epfield.μ[idx],epfield.s[idx],epfield.av[idx],epfield.va[idx], epfield, returnstatus, iter)
    end
end

function maxent_ep(metnet::MetNet; kwargs...)
    return maxent_ep(metnet.S, metnet.b, metnet.lb, metnet.ub; kwargs...)
end
