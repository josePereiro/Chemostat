## ------------------------------------------------------------------
# medium concentration
# s = c + u*ξ if u > 0 means outtake
function medium_conc(model::MetNet, out::AbstractOut, xi::Real, rxn::IDER_TYPE;
        intake_info = model.intake_info, 
        sense = 1, 
        force_pos = true,
        feed_c = () -> get(get(intake_info, rxn, Dict()), "c", 0.0)
    )
    rxn = rxns(model, rxn)
    u = sense * av(model, out, rxn)
    uerr = sense * sqrt(va(model, out, rxn))
    c = feed_c()
    s = c + u * xi
    serr = uerr * xi
    return force_pos ? (max(s, 0.0), serr) : (s, serr)
end

function medium_conc(bundle::ChstatBundle, xi::Real, k, rxn::IDER_TYPE; 
        kwargs...)
    model = bundle[xi, :net];
    out = bundle[xi, k]
    return medium_conc(model, out, xi, rxn; kwargs...)
end

function medium_conc(bundle::ChstatBundle, xi::Real, k1, k2, rxn::IDER_TYPE; 
    kwargs...)
    model = bundle[xi, :net];
    out = bundle[xi, k1, k2]
    return medium_conc(model, out, xi, rxn; kwargs...)
end
