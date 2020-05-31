struct HRout{T<:Real} <: AbstractOut
    av::Vector{T} # mean of the samples
    va::Vector{T} # variance of the samples
    hists::Vector{<:AbstractHistogram} # an histogram representation of the samples
    N::Integer # Number of flxs
    W::Integer # Number of samples
    hrsamples::AbstractArray # HR Samples (nsamples x fluxes)
end

function HRout(model::MetNet, hrsamples::AbstractArray; drop_samples = false, 
    step_factor = 0.01)
    isempty(hrsamples) && error("'hrsamples' is empty!!!")
    av = [mean(sample) for sample in eachcol(hrsamples)]
    va = [var(sample) for sample in eachcol(hrsamples)]
    hists = Histogram[]
    for (i, sample) in hrsamples |> eachcol |> enumerate
        step = (maximum(sample) - minimum(sample)) * step_factor
        if step != 0.0
            bins = model.lb[i]:step:model.ub[i]
            hist = fit(Histogram, sample, bins)
        else
            hist = fit(Histogram, sample)
        end
        push!(hists, hist)
    end
    W, N = size(hrsamples)
    hrsamples = drop_samples ? [] : hrsamples
    return HRout(av, va, hists, N, W, hrsamples)
end

HRout(hrout::HRout; drop_samples = false) = HRout(copy(hrout.av), copy(hrout.va), 
                deepcopy(hrout.hists), hrout.N, hrout.W, 
                drop_samples ? [] : copy(hrout.hrsamples))