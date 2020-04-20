const warn_color = :yellow
const info_color = :blue
const error_color = :red

"""
    Print useful info about the given metnet.
"""
function summary(metnet::MetNet)
    print("SUMMARY (color code: ")
    printstyled("warning", color = warn_color)
    printstyled(", info", color = info_color)
    printstyled(", error", color = error_color)
    println(")")
    printstyled("model size: $(size(metnet))", "\n", color = info_color)
    summary_rxn_bounds(metnet)
end

function summary_rxn_bounds(metnet::MetNet)
    length(metnet.lb) != metnet.N && (printstyled(
        "lb ($(length(metnet.lb))) != N ($(metnet.N)), dimention missmatch", 
        "\n", color = error_color); return)

    length(metnet.ub) != metnet.N && (printstyled(
        "ub ($(length(metnet.ub))) != N ($(metnet.N)), dimention missmatch", 
        "\n", color = error_color); return)

    all(iszero(metnet.lb)) && printstyled("all lb are zero", "\n", color = info_color)
    all(iszero(metnet.ub)) && printstyled("all ub are zero", "\n", color = error_color)
    flush(stdout);
    all(iszero(metnet.lb)) && all(iszero(metnet.ub)) && return
    printstyled("lb unique elments: ", sort!(unique(metnet.lb)), "\n",color = info_color)
    printstyled("ub unique elments: ", sort!(unique(metnet.ub)), "\n", color = info_color)

    for (i, rxn) in enumerate(metnet.rxns)
        lb = metnet.lb[i]
        ub = metnet.ub[i]
        lb > ub && printstyled("rxn($i): ($rxn), lb ($lb) > ub ($ub)", 
            "\n", color = error_color)
        lb > 0.0 && printstyled("rxn($i): ($rxn), lb ($lb) > 0.0", 
            "\n", color = warn_color);  
        lb == ub && printstyled("rxn($i): ($rxn), lb ($lb) == ub ($ub)", 
            "\n", color = info_color); 
        (isrev(metnet, i) ⊻ metnet.rev[i]) && printstyled("rxn($i): ($rxn), rev and bounds missmatch", 
            "\n", color = error_color)
        flush(stdout)
    end

    revscount(metnet) > 0 && printstyled("revscount: $(revscount(metnet))", "\n", color = error_color)
    fwdscount(metnet) > 0 && printstyled("fwds: $(fwdscount(metnet))", "\n", color = info_color)
    bkwdscount(metnet) > 0 && printstyled("bkwds: $(bkwdscount(metnet))", "\n", color = error_color)
    blockscount(metnet) > 0 && printstyled("blocks: $(blockscount(metnet))", "\n", color = error_color)
    fixxedscount(metnet) > 0 && printstyled("fixxed: $(fixxedscount(metnet))", "\n", color = info_color)

    return nothing
end