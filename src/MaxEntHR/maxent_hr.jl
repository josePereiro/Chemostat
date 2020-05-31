# TODO generalize this methods for maxent probability distriburions that depends not only in biomass

maxent_P(v::Real, β::Real, lb::Real, ub::Real) = β == 0 ? 1/(ub - lb) : β/(exp(β*(ub-v)) - exp(β*(lb - v)))

"""
    Resample over the politope but using the maxent distribution
    P ~ exp(beta * biomass)
    Returns the sub samples indexs

    - hrout: the solution of hr at beta cero
    - nsamples: number of samples to collect
    - maxiter: max number of iterations before return, even in nsamples is not collected
    - divs: number of divitions for discretizing biomass space to compute Z. (sorry if this is uggly)
"""
function maxent_hr_idxs(metnet::MetNet, hrout_β0::HRout, biomider, β; 
        nsamples = size(hrout_β0, 1),
        maxiter = 100 * nsamples,
        verbose = true,
        upfrec = 100
    )
    nsamples = floor(Int, nsamples)
    maxiter = floor(Int, maxiter)
    upfrec = floor(Int, nsamples/upfrec)
    verbose && (print("it 0: samples [0 / $(nsamples)]        \r"); flush(stdout))

    if isempty(hrout_β0.hrsamples)
        error("HRout hrsamples is empty, samples was probably droped!!!")
    end
    
    obj_idx = rxnindex(metnet, biomider)
    # The metnet must be preprocess (MetabolicEP.preprocess)
    # to ensure this is a the valid maximum (the same of FBA)
    obj_ub = metnet.ub[obj_idx] 
    obj_lb = metnet.lb[obj_idx]
    # TODO see TODO in the top
    # Only the biomass value is important
    obj_samples_β0 = view(hrout_β0.hrsamples, :, obj_idx)
    samples_β0_iter = 1:length(obj_samples_β0)

    maxent_samples_idxs = []
    for it in 1:maxiter
        rindx = rand(samples_β0_iter) # pick a random sample index
        robj_sample = obj_samples_β0[rindx] # random sample
        p_rsample = maxent_P(robj_sample, β, obj_lb, obj_ub) # maxent probability
        
        # Resampling
        if rand() < p_rsample
            push!(maxent_samples_idxs, rindx)
            verbose && it % upfrec == 0 && (print("it: $it  samples [$(length(maxent_samples_idxs)) / $(nsamples)]   \r"); flush(stdout))
        end

        if length(maxent_samples_idxs) == nsamples || it == maxiter
            verbose && (print("it: $it  samples [$(length(maxent_samples_idxs)) / $(nsamples)]   \r"); flush(stdout))
            break
        end
    end
    return maxent_samples_idxs
end

function maxent_hr(metnet::MetNet, hrout_β0::HRout, obj_ider::IDER_TYPE, β; 
        drop_samples = true, 
        nsamples = size(hrout_β0, 1),
        maxiter = 100 * nsamples
    )
    if isempty(hrout_β0.hrsamples)
        error("HRout hrsamples is empty, samples was probably dropped!!!")
    end
    idxs = maxent_hr_idxs(metnet, hrout_β0, obj_ider, β; 
        nsamples = nsamples, maxiter = maxiter)
    return HRout(metnet, view(hrout_β0.hrsamples, idxs,:), drop_samples = drop_samples)
end


function maxent_hr(metnet::MetNet, obj_ider::IDER_TYPE = "", β::Real = 0.0;
        nsamples = 10_000, maxiter = 1e10, 
        drop_samples = true,
        verbose = true,
        upfrec = 100
        )
    
    upfrec = floor(Int, nsamples/upfrec)
    hr = HRAlg(metnet)
    samples = zeros(nsamples, size(metnet, 2))
    
    verbose && (print("it 0: samples [0 / $(nsamples)]        \r"); flush(stdout))
    if β == 0
        for (it, sample) in enumerate(eachrow(samples))
            sample .= hrsample!(hr)
            verbose && it % upfrec == 0 && (print("it: $it  samples [$it / $(nsamples)]   \r"); flush(stdout))
        end
    else
    
        obj_idx = rxnindex(metnet, obj_ider)
        obj_ub = metnet.ub[obj_idx] 
        obj_lb = metnet.lb[obj_idx]
        
        spos = 1
        for it in 1:floor(Int, maxiter)
            sample = hrsample!(hr)
            obj_val = sample[obj_idx]
            p_sample = maxent_P(obj_val, β, obj_lb, obj_ub) # maxent probability
            if rand() < p_sample
                verbose && spos % upfrec == 0 && 
                    (print("it: $it  samples [$(spos -1) / $(nsamples)]   obj_val: $(obj_val)         \r"); flush(stdout))
                samples[spos, :] .= sample
                spos += 1
            end
            spos > nsamples && break
        end
    end
    return HRout(metnet, samples, drop_samples = drop_samples)
end