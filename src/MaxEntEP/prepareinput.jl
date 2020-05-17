# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# TODO remove Y unused parameter
# TODO remove alpha Inf code
function prepareinput(K, Y, lb, ub, beta, verbose, solution, expval)

    M,N = size(K)
    M < N || @warn("M = $M ≥ N = $N")
    all(lb .<= ub) || error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

    verbose && println(stderr, "Analyzing a $M × $N stoichiometric matrix.")

    scalefact = zero(eltype(K))

    updatefunction = if beta == Inf
        eponesweepT0!
    else
        eponesweep!
    end

    scalefact = max(maximum(abs.(lb)), maximum(abs.(ub)))
    if solution === nothing
        epfield = EPFields(N,expval,scalefact)
    else
        epfield = deepcopy(solution.sol) # preserve the original solution!
    end

    return updatefunction, scalefact, epfield
end