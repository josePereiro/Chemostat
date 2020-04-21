av(epout) = epout.av
av(metnet::MetNet, epout::EPout, ider) = epout.av[rxnindex(metnet, ider)]
av(epout::EPout, metnet::MetNet, ider) = av(metnet, epout, ider)

va(epout) = epout.va
va(metnet::MetNet, epout::EPout, ider) = epout.va[rxnindex(metnet, ider)]
va(epout::EPout, metnet::MetNet, ider) = va(metnet, epout, ider)

μ(epout) = epout.av
μ(metnet::MetNet, epout::EPout, ider) = epout.μ[rxnindex(metnet, ider)]
μ(epout::EPout, metnet::MetNet, ider) = μ(metnet, epout, ider)

σ(epout) = epout.av
σ(metnet::MetNet, epout::EPout, ider) = epout.σ[rxnindex(metnet, ider)]
σ(epout::EPout, metnet::MetNet, ider) = σ(metnet, epout, ider)