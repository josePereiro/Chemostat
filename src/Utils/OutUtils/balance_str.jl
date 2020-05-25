function balance_str(model::MetNet, out::AbstractOut, ider::IDER_TYPE; digits = 50)
    meti = metindex(model, ider)
    rxns = model.rxns[met_rxns(model, ider)] |> sort
    b_str = []
    b = 0.0
    for rxn in rxns
        s = round(S(model, meti, rxn), digits = digits)
        f = round(av(model, out, rxn), digits = digits)
        sf = round(s*f, digits = digits)
        b += s*f
        push!(b_str, "($s*$f=$sf)$(rxn)")
    end
    return "$(model.mets[meti]): " * join(b_str, " + ") * " [tot: $(round(b, digits = digits))]" * " == " * 
                "$(round(model.b[meti], digits = digits))"
end
balance_str(model, out::AbstractOut, iders::Vector; digits = 50) = 
    [balance_str(model, out, ider; digits = digits) for ider in iders]