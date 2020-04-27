# TODO generalize this methods for maxent probability distriburions that depends not only in biomass

maxent_fun(biom::Real, biom_ub::Real, β::Real)::Real = exp(-β*(biom_ub - biom))

function maxent_biomass_Z(biom_lb::Real, biom_ub::Real, β::Real; divs::Integer = 1_000_000)::Real
    Z = 0.0
    biom_step = (biom_ub - biom_lb)/divs
    for biom in biom_lb:biom_step:biom_ub
        Z += maxent_fun(biom, biom_ub, β)*biom_step
    end
    return Z
end

"""
    Resample over the politope but using the maxent distribution
    P ~ exp(beta * biomass)
    Returns the sub samples indexs

    - nsamples: number of samples to collect
    - max_iters: max number of iterations before return, even in nsamples is not collected
    - divs: number of divitions for discretizing biomass space to compute Z. (sorry if this is uggly)
"""
function maxent_hrsamples_idxs(metnet::MetNet, hrout::HRout, biomider, β; 
        nsamples = size(hrout, 1),
        max_iters = 100 * nsamples,
        divs = 1_000_000 # for computing Z
    )

    if isempty(hrout.hrsamples)
        error("HRout hrsamples is empty, samples was probably droped!!!")
    end
    
    W, N = size(hrout)
    biom_idx = rxnindex(metnet, biomider)
    # The metnet must be preprocess (MetabolicEP.preprocess)
    # to ensure this is a the valid maximum (the same of FBA)
    biom_ub = metnet.ub[biom_idx] 
    biom_lb = metnet.lb[biom_idx]
    # TODO see TODO in the top
    # Only the biomass value is important
    biom_samples = hrout.hrsamples[:, biom_idx]
    
    Z = maxent_biomass_Z(biom_lb, biom_ub, β, divs = divs)
    
    maxent_samples_idxs = []
    for it in 1:max_iters
        rindx = rand(1:W) # pick a random sample index
        rbiom_sample = biom_samples[rindx] # random sample
        p_rsample = maxent_fun(rbiom_sample, biom_ub, β)/Z # maxent probability
        
        # Resampling
        rand() < p_rsample && push!(maxent_samples_idxs, rindx)

        length(maxent_samples_idxs) == nsamples && break
    end
    return maxent_samples_idxs
end

function maxent_hrsamples(metnet, hrout::HRout, biomider, β; 
        drop_samples = true, 
        nsamples = size(hrout, 1),
        max_iters = 100 * nsamples,
        divs = 1_000_000
    )
    if isempty(hrout.hrsamples)
        error("HRout hrsamples is empty, samples was probably dropped!!!")
    end
    idxs = maxent_hrsamples_idxs(metnet, hrout, biomider, β; 
        nsamples = nsamples, max_iters = max_iters, divs = divs)
    return HRout(hrout.hrsamples[idxs,:], drop_samples)
end