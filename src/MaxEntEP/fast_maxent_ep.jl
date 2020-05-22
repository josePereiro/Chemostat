"""
    Make MaxEnt-EP but applaying MaxEnt starting from a EP solution at beta zero.
    It is important to pass a epout solution computed with beta zero.
    It is faster than the fully integrated version 'maxent_ep' but least accurate.
    It is used for secundary computations as estimating a beta which lead to a 
    given objective value and so use this beta for the fully integrated version of the algorithm
"""
function fast_maxent_ep(model::MetNet, epout_β0::EPout, α::Real,  βv::Vector)

    #= Anna's code scales fluxes to [-1,1], but it 
    doesn't rescale a,d back to the original range =#
    maxflux = max(maximum(abs.(model.lb)), maximum(abs.(model.ub)))
    a = epout_β0.sol.a * maxflux;
    d = epout_β0.sol.b * maxflux^2;
    #= a, d are the mean and variance of the univariate Gaussians
    used as site approximations in the EP =#

    #= v, Σ are the mean vector and covariance matrix of Q,
    the full multivariate Gaussian (in Braunstein et al paper) =#
    S = model.S
    b = model.b
    D = Diagonal(1.0 ./ d)
    invΣ = Matrix(α*S'*S + D)
    Σ = inv(invΣ)
    v = Σ*(α*S'*b + D*a)

    #= From Cossios papar (see README) =#
    w = v + Σ * βv
    # wn = (d .* w - diag(Σ) .* a) ./ (d .- diag(Σ)) # mean of the non-truncated marginals
    # Σnn = d .* diag(Σ) ./ (d .- diag(Σ)) # variances of the non-truncated marginals

    tns = Truncated.(Normal.(wn, sqrt.(Σnn)), model.lb, model.ub)

    return EPout(wn, Σnn, mean.(tns), var.(tns), nothing, epout_β0.status)
end