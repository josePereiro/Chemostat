function EPModel(S::AbstractArray{T,2}, b::Array{T,1}, lb::Array{T,1}, ub::Array{T,1};
    alpha::Real=1e7,      # inverse temperature
    beta_vec::Vector{T} = T[], # maxent inverse temperature vector
    solution::Union{EPout{T}, Nothing} = nothing,  # start from a solution
    expval=nothing # fix posterior probability experimental values for std and mean
) where {T<:Real}

# Some checks
M,N = size(S)
M > N && @warn("M = $M ≥ N = $N")
all(lb .<= ub) || error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

# The scalefactor is just the maximum absolute bound (lb or ub).
scalefact = max(maximum(abs.(lb)), maximum(abs.(ub)))

# Create EPField, If a solution is not given, the EPfield will be fresh
epfield = isnothing(solution) ? EPFields(N, expval, eltype(S)) : 
    deepcopy(solution.sol) # preserve the original solution!

# making a local copy to rescale
llb = copy(lb) 
lub = copy(ub)
lb = copy(b)

#=
Scale down μ, s, av, va of epfield and lub, llb and Y using the 
previous computed scalefactor.
If epfield is fresh, it only will have effect on lub, llb and Y
=#
scaleepfield!(epfield, 1.0/scalefact, lub, llb, lb) # scaling fields in [0,1]

epmat = alpha < Inf ? EPMat(S, lb, llb, lub, alpha) : EPMatT0(S, lb, llb, lub)

# One iteration of EP
updatealg! = alpha == Inf ? eponesweepT0! : eponesweep!

beta_vec = prepare_βv(epmat, beta_vec)

return EPModel(scalefact, updatealg!, epfield, epmat, alpha, beta_vec)

end

EPModel(metnet::MetNet; kwargs...) = EPModel(metnet.S, metnet.b, metnet.lb, metnet.ub; kwargs...)