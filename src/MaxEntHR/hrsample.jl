# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function hrsample!(hr::HRAlg)
    # pick a random direction
    hr.v .= hr.base * randn(hr.k)
    dx = hr.v/norm(hr.v)
    # compute intersection
    l = maximum(min((hr.lb[i]-hr.x[i])/dx[i], (hr.ub[i]-hr.x[i])/dx[i]) for i=1:hr.n if dx[i] != 0)
    u = minimum(max((hr.lb[i]-hr.x[i])/dx[i], (hr.ub[i]-hr.x[i])/dx[i]) for i=1:hr.n if dx[i] != 0)
    # find a random point in the intersection
    t = l+(u-l)*rand()
    hr.x .+= t * dx
    return hr.x
end

function hrsample!(hr::HRAlg, nsamples::Integer)
    samples = zeros(nsamples, hr.n)
    for sample in eachrow(samples)
        sample .= hrsample!(hr)
    end
    return samples
end