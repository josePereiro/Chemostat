function metabolicEP(metnet::MetNet; include_epmat = false, kwargs...)
    epout, epmat = metabolicEP(metnet.S, metnet.b, metnet.lb, metnet.ub; kwargs...)
    return include_epmat ? (epout, epmat) : epout;
end

#TODO make pull request about inconcistence with name K -> S and Y -> b

"""
This is an overwitten version of the method taken from:
Alfredo Braunstein, Anna Paola Muntoni, and Andrea Pagnani, “An Analytic Approximation of the Feasible 
Space of Metabolic Networks,” Nature Communications 8 (April 6, 2017): https://doi.org/10.1038/ncomms14915.
Package metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

The only change is that now the method returns also the epmat object, required for
computing the MaxEnt distribution

Original doc:

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
- `beta` (inverse temperature::``Real``): default 10^7;  the zero temperature algorithm is run setting ``beta=Inf``.
- `verbose` (``true`` or ``false``): default ``true``
- `damp` (∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield): default 0.9
- `epsconv` (convergence criterion): default 1e-6
- `maxiter` (maximum number of iterations): default 2000
- `maxvar`  (threshold on maximum variance): default 1e50
- `minvar`  (threshold on minimum variance): default 1e-50
- `solution` (start from solution. Is of type ``EPout``): default: ``nothing``
- `expval` (fix to posterior probability of mean and/or variance to values): default ``nothing``. expval can be either at ``Tuple{Float64,Float64,Int}`` or a ``Vector{Tuple{Float64,Float64,Int}}``. Values can be fixed as``expval=(0.2,0.4,4)`` meaning that for flux index 4 the mean is set to 0.2 and the variance to 0.4. Fixing more values ``expval=[(0.2, 0.3, 4), (0.4, nothing, 5)]``: in this case, we fix the posterior of flux 4 to 0.2 (mean) and 0.3 (variance), while for flux 5 we fix the mean to 0.4 and we keep the variance free.

"""
function metabolicEP(S::AbstractArray{T,2}, b::Array{T,1}, lb::Array{T,1}, ub::Array{T,1};
                     beta::Real=1e7,      # inverse temperature
                     verbose::Bool=true,  # output verbosity
                     damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
                     epsconv::Real=1e-6,  # convergence criterion
                     maxiter::Int=2000,   # maximum iteration count
                     maxvar::Real=1e50,   # maximum numerical variance
                     minvar::Real=1e-50,  # minimum numerical variance
                     solution::Union{MetabolicEP.EPout{T},Nothing}=nothing,  # start from a solution
                     expval=nothing # fix posterior probability experimental values for std and mean
                     ) where T<:Real

    llb = copy(lb) # making  a local copy to rescale
    lub = copy(ub)
    updatealg,scalefact,epfield = MetabolicEP.prepareinput(S,b,llb,lub,beta,verbose,solution,expval)

    MetabolicEP.scaleepfield!(epfield,lub,llb,b,1.0/scalefact) # rescaling fields in [0,1]

    epalg = MetabolicEP.EPAlg(beta, minvar, maxvar, epsconv, damp, maxiter,verbose)
    epmat = beta < Inf ? MetabolicEP.EPMat(S,b,llb, lub, beta) : MetabolicEP.EPMatT0(S,b,llb, lub)
    returnstatus = MetabolicEP.epconverge!(epfield,epmat,epalg, updatealg)
    MetabolicEP.scaleepfield!(epfield,lub,llb,b,scalefact)
    if beta < Inf
        return  MetabolicEP.EPout(epfield.μ,epfield.s, epfield.av, epfield.va, epfield, returnstatus), epmat
        # in Anna's original code epmat is not returned
    else
        idx = epmat.idx
        return  MetabolicEP.EPout(epfield.μ[idx],epfield.s[idx],epfield.av[idx],epfield.va[idx], epfield, returnstatus), epmat
        # in Anna's original code epmat is not returned
    end
end