function EPModel(S::AbstractArray{T,2}, b::AbstractArray{T}, lb::AbstractArray{T}, ub::AbstractArray{T};
        alpha::Real=1e7,      # inverse temperature
        beta_vec::AbstractVector{T} = spzeros(T, size(S, 2)), # maxent inverse temperature vector
        solution::Union{EPout{T}, Nothing} = nothing,  # start from a solution
        expval = nothing # fix posterior probability experimental values for std and mean
    ) where {T<:Real}

    # Some checks
    M, N = size(S)
    M > N && @warn("M = $M ≥ N = $N")
    any(lb .> ub) && error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

    # The scalefactor is just the maximum absolute bound (lb or ub).
    scalefact = max(maximum(abs.(lb)), maximum(abs.(ub)))

    # Create EPFields, If a solution is not given, the EPfields will be fresh
    epfields = isnothing(solution) ? EPFields(N, expval, eltype(S)) : 
        deepcopy(solution.sol) # preserve the original solution!

    # making a local copy to rescale
    lb, ub, b = copy.([lb, ub, b]) 

    #=
    Scale down μ, s, av, va of epfields and ub, lb and Y using the 
    previous computed scalefactor.
    If epfields is fresh, it only will have effect on ub, lb and Y
    =#
    scaleepfield!(inv(scalefact), epfields, ub, lb, b) # scaling fields in [0,1]

    epmat = alpha < Inf ? EPMat(S, b, lb, ub, alpha) : EPMatT0(S, b, lb, ub)

    # One iteration of EP
    updatealg! = alpha == Inf ? eponesweepT0! : eponesweep!

    beta_vec = prepare_βv(N, beta_vec)

    return EPModel{T}(scalefact, updatealg!, epfields, epmat, alpha, beta_vec, Dict())

end

EPModel(metnet::MetNet; kwargs...) = EPModel(metnet.S, metnet.b, metnet.lb, metnet.ub; kwargs...)