function balance_str(model::MetNet, out, ider::IDER_TYPE; digits = 50)
    meti = metindex(model, ider)
    rxns = met_rxns(model, ider)
    b_str = []
    b = 0.0
    for rxn in rxns
        s = model.S[meti, rxn]
        f = av(model, out, rxn)
        b += s*f
        push!(b_str, "$(model.rxns[rxn])($s*$(round(f, digits = digits)) = $(round(s*f, digits = digits)))")
    end
    return "$(model.mets[meti]): " * join(b_str, " + ") * " == " * 
                "$(round(model.b[meti], digits = digits))" * " [$(round(b, digits = digits))]"
end
balance_str(model, out, iders::Vector; digits = 50) = 
    [balance_str(model, out, ider; digits = digits) for ider in iders]