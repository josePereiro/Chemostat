# From Cossios MaxEntChemostat2018
function maxent_metabolicEP(model::MetNet, epout::EPout, epmat::AbstractEPMat, obj_ider, β0; verbose = true)
    
    # TODO see how to implement this with beta = Inf
    epmat isa EPMatT0 && error("Not implemented for metabolicEP, beta = Inf")
    
    #= Anna's code scales fluxes to [-1,1], but it 
    doesn't rescale a,d back to the original range =#
    maxflux = max(maximum(abs.(model.lb)), maximum(abs.(model.ub)))
    a = epout.sol.a * maxflux;
    d = epout.sol.b * maxflux^2;
    #= a, d are the mean and variance of the univariate Gaussians
    used as site approximations in the EP =#
    
    #= Similarly, Anna's code never rescales the 
    contents of epmat =#
    v = epmat.v * maxflux
    Σ = epmat.invKKPD * maxflux^2
    iΣ = epmat.KKPD / maxflux^2
    #= v, Σ are the mean vector and covariance matrix of Q,
    the full multivariate Gaussian (in Braunstein et al paper) =#

    #= The mean and variance of the univariate marginals of Q(n)    
    for each flux. Note that this does not include truncation.
    These guys are re-scaled in Anna's code, so no need to do that
    here again. =#
    μ = epout.μ
    ν = epout.σ

    # I dont want to introduce a new intermediate data type
    # I'll try to return a new EPout
    
    # 
    obj_ider = rxnindex(model, obj_ider)
    α = spzeros(size(Σ, 1))
    α[obj_ider] = β0

    w = v + Σ * α
    Σnn = d .* diag(Σ) ./ (d .- diag(Σ)) # variances of the non-truncated marginals
    wn = (d .* w - diag(Σ) .* a) ./ (d .- diag(Σ)) # mean of the non-truncated marginals
    tns = Truncated.(Normal.(wn, sqrt.(Σnn)), model.lb, model.ub)
    
    # TODO check deeper the implecations of passing the same 'sol'
    return typeof(epout)(wn, Σnn, mean.(tns), var.(tns), epout.sol, epout.status)
    
end