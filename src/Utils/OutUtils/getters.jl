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
hists(model::MetNet, hrout::HRout, ider::IDER_TYPE) = 
    (ider = rxnindex(model, ider); hrout.hists[ider])

hrsamples(hrout::HRout) = hrout.hrsamples
hrsamples(model::MetNet, hrout::HRout, ider::IDER_TYPE) = 
    (ider = rxnindex(model, ider); hrout.hrsamples[:, ider])

# Commons getter interface
for fun in [av, va, μ, σ]
    fun_name = string(nameof(fun))

    eval(Meta.parse(
        """$(fun_name)(metnet::MetNet, out, ider::IDER_TYPE) = 
                $(fun_name)(out)[rxnindex(metnet, ider)]"""))
    eval(Meta.parse(
        """$(fun_name)(metnet::MetNet, out, iders::Vector) = 
                [$(fun_name)(metnet, out, ider) for ider in iders]"""))
    eval(Meta.parse(
        """$(fun_name)(metnet::MetNet, outs::Vector, ider::IDER_TYPE) =
                [$(fun_name)(metnet, out, ider) for out in outs]"""))
    eval(Meta.parse(
        """$(fun_name)(metnets::Vector, outs::Vector, ider::IDER_TYPE) = 
                [$(fun_name)(metnet, out, ider) for (metnet, out) in zip(metnets, outs)]"""))
end