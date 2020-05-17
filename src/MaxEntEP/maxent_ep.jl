# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function maxent_ep(K::AbstractArray{T,2}, Y::Array{T,1}, lb::Array{T,1}, ub::Array{T,1};
    α::Real=1e7,      # inverse temperature
    βv::Vector{T} = T[], # maxent inverse temperature vector
    verbose::Bool=true,  # output verbosity
    damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
    epsconv::Real=1e-6,  # convergence criterion
    maxiter::Int=2000,   # maximum iteration count
    maxvar::Real=1e50,   # maximum numerical variance
    minvar::Real=1e-50,  # minimum numerical variance
    solution::Union{EPout{T}, Nothing}=nothing,  # start from a solution
    expval=nothing # fix posterior probability experimental values for std and mean
    ) where T<:Real


    llb = copy(lb) # making  a local copy to rescale
    lub = copy(ub)

    #= 
    Create EPField, set scale factor and return updatealg (updatefunction),
    depending on beta (eponesweepT0! or eponesweep).
    The scalefactor is just the maximum absolute bound (lb or ub).
    If a solution is not given, the EPfield will be fresh
    =#
    updatealg,scalefact,epfield = prepareinput(K,Y,llb,lub, α, verbose,solution,expval)

    #=
    Scale down μ, s, av, va of epfield and lub, llb and Y using the 
    previous computed scalefactor.
    If epfield is fresh, it only will have effect on lub, llb and Y
    =#
    scaleepfield!(epfield,lub,llb,Y,1.0/scalefact) # rescaling fields in [0,1]

    epmat = α < Inf ? EPMat(K, Y, llb, lub, α) : EPMatT0(K,Y,llb, lub)
    βv = prepare_βv(epmat, βv)
    epalg = EPAlg(α, βv, minvar, maxvar, epsconv, damp, maxiter, verbose)

    #= Here is were all the work is done, this function will 
    call updatealg till convergence or maxiter is reached =#
    returnstatus = epconverge!(epfield, epmat, epalg, updatealg)

    #= Scale back μ, s, av, va of epfield and lub, llb and Y =#
    scaleepfield!(epfield,lub,llb,Y,scalefact)
    if α < Inf
        return  EPout(epfield.μ,epfield.s, epfield.av, epfield.va, epfield, returnstatus)
    else
        idx = epmat.idx
        return  EPout(epfield.μ[idx],epfield.s[idx],epfield.av[idx],epfield.va[idx], epfield, returnstatus)
    end
end