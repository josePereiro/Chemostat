struct HRout{T<:AbstractFloat}
    av::Vector{T} # mean of the samples
    va::Vector{T} # variance of the samples
    hists::Vector{<:AbstractHistogram} # an histogram representation of the samples
    N::Integer # Number of flxs
    W::Integer # Number of samples
    hrsamples::AbstractArray # HR Samples (nsamples x fluxes)
end

function HRout(hrsamples::AbstractArray, drop_samples = false)
    isempty(hrsamples) && error("'hrsamples' is empty!!!")
    av = [mean(sample) for sample in eachcol(hrsamples)]
    va = [var(sample) for sample in eachcol(hrsamples)]
    hists = [fit(Histogram, sample, nbins = 100) for sample in eachcol(hrsamples)]
    W, N = size(hrsamples)
    hrsamples = drop_samples ? [] : hrsamples
    return HRout(av, va, hists, N, W, hrsamples)
end