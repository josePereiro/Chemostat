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

# Commons
av(metnet::MetNet, out, ider) = av(out)[rxnindex(metnet, ider)]
av(metnet::MetNet, outs::Vector, ider) =
    [av(metnet, out, ider) for out in outs]
av(metnets::Vector, outs::Vector, ider) = 
    [av(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

va(metnet::MetNet, out, ider) = va(out)[rxnindex(metnet, ider)]
va(metnet::MetNet, outs::Vector, ider) = 
    [va(metnet, out, ider) for out in outs]
va(metnets::Vector, outs::Vector, ider) = 
    [va(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

μ(metnet::MetNet, out, ider) = μ(out)[rxnindex(metnet, ider)]
μ(metnet::MetNet, outs::Vector, ider) = 
    [μ(metnet, out, ider) for out in outs]
μ(metnets::Vector, outs::Vector, ider) = 
    [μ(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

σ(metnet::MetNet, out, ider) = σ(out)[rxnindex(metnet, ider)]
σ(metnet::MetNet, outs::Vector, ider) = 
    [σ(metnet, out, ider) for out in outs]
σ(metnets::Vector, outs::Vector, ider) = 
    [σ(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]

