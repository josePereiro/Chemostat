function balance_str(model, ider; digits = 50)
    meti = metindex(model, ider)
    rxns = model.rxns[met_rxns(model, ider)] |> sort
    b_str = []
    for rxn in rxns
        s = round(S(model, meti, rxn), digits = digits)
        push!(b_str, "($s)$rxn")
    end
    return "$(model.mets[meti]): " * join(b_str, " + ") * " == " * 
                "$(round(model.b[meti], digits = digits))"
end