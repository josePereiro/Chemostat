"""
    Print useful info about the given metnet.
"""
function summary(metnet::MetNet)
    summary_rxn_bounds(metnet)
end

function summary_rxn_bounds(metnet::MetNet)
    println("Bounds")

    length(metnet.lb) != metnet.N && (printstyled(
        "lb ($(length(metnet.lb))) != N ($(metnet.N)), dimention missmatch", 
        "\n", color = :red); return)

    length(metnet.ub) != metnet.N && (printstyled(
        "ub ($(length(metnet.ub))) != N ($(metnet.N)), dimention missmatch", 
        "\n", color = :red); return)

    all(iszero(metnet.lb)) && printstyled("all lb are zero", "\n", color = :blue)
    all(iszero(metnet.ub)) && printstyled("all ub are zero", "\n", color = :blue)
    flush(stdout);
    all(iszero(metnet.lb)) && all(iszero(metnet.ub)) && return
    printstyled("lb unique elments: ", sort!(unique(metnet.lb)), "\n",color = :blue)
    printstyled("ub unique elments: ", sort!(unique(metnet.ub)), "\n", color = :blue)

    for (i, rxn) in enumerate(metnet.rxns)
        lb = metnet.lb[i]
        ub = metnet.ub[i]
        lb > ub && printstyled("rxn($i): ($rxn), lb ($lb) > ub ($ub)", 
            "\n", color = :red)
        lb > 0.0 && printstyled("rxn($i): ($rxn), lb ($lb) > 0.0", 
            "\n", color = :blue);  
        lb == ub && printstyled("rxn($i): ($rxn), lb ($lb) == ub ($ub)", 
            "\n", color = :blue); 
        ((lb < 0.0 && ub > 0.0) ⊻ metnet.rev[i]) && printstyled("rxn($i): ($rxn), rev and bounds missmatch", 
            "\n", color = :red)
        flush(stdout)
    end
end