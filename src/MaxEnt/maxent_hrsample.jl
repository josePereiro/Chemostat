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
"""
function maxent_hrsamples_idxs(metnet, hrsamples, β, biomider; divs = 1_000_000)
    biom_idx = rxnindex(metnet, biomider)
    
    biom_ub = metnet.ub[biom_idx]
    biom_lb = metnet.lb[biom_idx]
    Z = maxent_biomass_Z(biom_lb, biom_ub, β, divs = divs)
    
    maxent_samples_idxs = []
    for (i, sample) in enumerate(eachrow(hrsamples))
        biom = sample[biom_idx]
        p_biom = maxent_fun(biom, biom_ub, β)/Z
        
        # Resampling
        rand() < p_biom && push!(maxent_samples_idxs, i)
    end
    return maxent_samples_idxs
end

function maxent_hrsamples(metnet, biomider, β; 
    niter::Integer = 1, nsamples_per_iter::Integer = 1_000_000)
    M, N = size(metnet)
    # get a bigger matrix than necessary
    maxent_samples = zeros(nsamples_per_iter * niter, N) 
    next_sample_pos = 1;

    # TODO parallelize this loop, then set the default value of niter to something bigger
    for it in 1:niter
        hrsamples = hrsample(metnet; nsamples = nsamples_per_iter)
        samples_idxs = maxent_hrsamples_idxs(metnet, hrsamples, β, biomider;)

        curr_range = next_sample_pos:(next_sample_pos + length(samples_idxs) - 1)
        
        maxent_samples[curr_range,:] = hrsamples[samples_idxs, :]
        next_sample_pos += length(samples_idxs)
    end
    return maxent_samples[1:(next_sample_pos - 1), :]
end