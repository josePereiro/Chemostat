# TODO generalize this methods for maxent probability distriburions that depends not only in biomass

maxent_Z(β::Real, lb::Real, ub::Real) = β == 0 ? 1/(ub - lb) : ((exp(β*ub) - exp(β*lb))/β)
unsafe_maxent_P(x::Real, β::Real, Z::Real) = exp(β*x)/Z
maxent_P(x::Real, β::Real, lb::Real, ub::Real, Z = maxent_Z(β, lb, ub)) = 
    lb <= x <= ub ? unsafe_maxent_P(x, β, Z) : 0.0

"""
    Resample over the politope but using the maxent distribution
    P ~ exp(beta * biomass)
    Returns the sub samples indexs

    - nsamples: number of samples to collect
    - max_iters: max number of iterations before return, even in nsamples is not collected
    - divs: number of divitions for discretizing biomass space to compute Z. (sorry if this is uggly)
"""
# function maxent_hr_idxs(metnet::MetNet, hrout::HRout, biomider, β; 
#         nsamples = size(hrout, 1),
#         max_iters = 100 * nsamples
#     )

#     if isempty(hrout.hrsamples)
#         error("HRout hrsamples is empty, samples was probably droped!!!")
#     end
    
#     W, N = size(hrout)
#     biom_idx = rxnindex(metnet, biomider)
#     # The metnet must be preprocess (MetabolicEP.preprocess)
#     # to ensure this is a the valid maximum (the same of FBA)
#     biom_ub = metnet.ub[biom_idx] 
#     biom_lb = metnet.lb[biom_idx]
#     # TODO see TODO in the top
#     # Only the biomass value is important
#     biom_samples = hrout.hrsamples[:, biom_idx]
    
#     Z = maxent_Z(β, biom_lb, biom_ub)
#     !isfinite(Z) && error("Z Inf, maybe beta is too large or the bounds "*
#         "of the objective reaction are too wide. Consider use Utils.fva_preprocess")
    
#     maxent_samples_idxs = []
#     for it in 1:max_iters
#         rindx = rand(1:W) # pick a random sample index
#         rbiom_sample = biom_samples[rindx] # random sample
#         p_rsample = unsafe_maxent_P(rbiom_sample, β, Z) # maxent probability
        
#         # Resampling
#         rand() < p_rsample && push!(maxent_samples_idxs, rindx)

#         length(maxent_samples_idxs) == nsamples && break
#     end
#     return maxent_samples_idxs
# end

# function maxent_hr(metnet::MetNet, hrout::HRout, biomider, β; 
#         drop_samples = true, 
#         nsamples = size(hrout, 1),
#         max_iters = 100 * nsamples
#     )
#     if isempty(hrout.hrsamples)
#         error("HRout hrsamples is empty, samples was probably dropped!!!")
#     end
#     idxs = maxent_hr_idxs(metnet, hrout, biomider, β; 
#         nsamples = nsamples, max_iters = max_iters)
#     return HRout(metnet, hrout.hrsamples[idxs,:], drop_samples)
# end


function maxent_hr(metnet::MetNet, obj_ider::IDER_TYPE, β::Real;
        nsamples = 10_000, maxiter = 1e10, 
        drop_samples = true)
    
    hr = HRAlg(metnet)
    obj_idx = rxnindex(metnet, obj_ider)
    
    obj_ub = metnet.ub[obj_idx] 
    obj_lb = metnet.lb[obj_idx]

    # Only the biomass value is important
    Z = maxent_Z(β, obj_lb, obj_ub)
    !isfinite(Z) && error("Z Inf, maybe beta is too large or the bounds "*
        "of the objective reaction are too wide. Consider use Utils.fva_preprocess")
    
    samples = zeros(nsamples, size(metnet, 2))
    spos = 1
    for it in 1:maxiter
        sample = hrsample!(hr)
        obj_val = sample[obj_idx]
        p_sample = unsafe_maxent_P(obj_val, β, Z) # maxent probability
        rand() < p_sample && (samples[spos, :] .= sample, spos += 1)
        spos > nsamples && break
    end
    
    return HRout(metnet, samples, drop_samples)
end