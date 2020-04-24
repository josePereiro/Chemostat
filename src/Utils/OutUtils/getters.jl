# FBAout
av(fbaout::FBAout) = fbaout.v
va(fbaout::FBAout) = zeros(length(fbaout.v))
μ(fbaout::FBAout) = fbaout.v
σ(fbaout::FBAout) = zeros(length(fbaout.v))

# EPout
av(epout::EPout) = epout.av
va(epout::EPout) = epout.va
μ(epout::EPout) = epout.μ
σ(epout::EPout) = epout.σ

# HRout
av(hrout::HRout) = hrout.av
va(hrout::HRout) = hrout.va
μ(hrout::HRout) = hrout.av
σ(hrout::HRout) = hrout.va

hists(hrout::HRout) = hrout.hist
hists(model::MetNet, hrout::HRout, ider) = 
    (ider = rxnindex(model, ider); hrout.hists[ider])

hrsamples(hrout::HRout) = hrout.hrsamples
hrsamples(model::MetNet, hrout::HRout, ider) = 
    (ider = rxnindex(model, ider); hrout.hrsamples[:, ider])

# Boundle
va_ep(boundle::ChstatBoundle, ξ, β, ider) = 
    va(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)
av_ep(boundle::ChstatBoundle, ξ, β, ider) = 
    av(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)
μ_ep(boundle::ChstatBoundle, ξ, β, ider) = 
    μ(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)
σ_ep(boundle::ChstatBoundle, ξ, β, ider) = 
    σ(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)

va_hr(boundle::ChstatBoundle, ξ, β, ider) = 
    va(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)
av_hr(boundle::ChstatBoundle, ξ, β, ider) = 
    av(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)
μ_hr(boundle::ChstatBoundle, ξ, β, ider) = 
    μ(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)
σ_hr(boundle::ChstatBoundle, ξ, β, ider) = 
    σ(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)

va_fba(boundle::ChstatBoundle, ξ, ider) = 
    va(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)
av_fba(boundle::ChstatBoundle, ξ, ider) = 
    av(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)
μ_fba(boundle::ChstatBoundle, ξ, ider) = 
    μ(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)
σ_fba(boundle::ChstatBoundle, ξ, ider) = 
    σ(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)



# Commons
av(metnet::MetNet, out, ider) = av(out)[rxnindex(metnet, ider)]
av(metnets::Vector, outs::Vector, ider) = 
    [av(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

va(metnet::MetNet, out, ider) = va(out)[rxnindex(metnet, ider)]
va(metnets::Vector, outs::Vector, ider) = 
    [va(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

μ(metnet::MetNet, out, ider) = μ(out)[rxnindex(metnet, ider)]
μ(metnets::Vector, outs::Vector, ider) = 
    [μ(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

σ(metnet::MetNet, out, ider) = σ(out)[rxnindex(metnet, ider)]
σ(metnets::Vector, outs::Vector, ider) = 
    [σ(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

