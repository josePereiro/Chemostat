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
    PRINT_MAX = 50 # TODO make this global?
    M, N = size(metnet)

    all(iszero(metnet.lb)) && all(iszero(metnet.ub)) && 
        (printstyled("lb and ub boths has only zero elements", "\n", color = error_color); return)
    
    # single checks
    for (id, val) in zip(["lb", "ub"], [metnet.lb, metnet.ub])
    
        length(val) != N && (printstyled(
            "$id ($(length(val))) != N ($N), dimention missmatch", 
            "\n", color = error_color); return)

        # length(metnet.ub) != N && (printstyled(
        #     "ub ($(length(metnet.ub))) != N ($N), dimention missmatch", 
        #     "\n", color = error_color); return)

        all(iszero(val)) && printstyled("all $id are zero", "\n", color = info_color)
        # all(iszero(metnet.ub)) && printstyled("all ub are zero", "\n", color = error_color)
    
        unique_ = sort!(unique(val))
        length(unique_) < PRINT_MAX ?
            printstyled("   $id $(length(unique_)) unique elment(s): ", unique_, "\n", color = info_color) :
            printstyled("   $id $(length(unique_)) unique elment(s): min: ", 
                first(unique_), " max: ", last(unique_), "\n", color = info_color)

    end

    # Counple checks
    line_count = 0
    for (i, rxn) in enumerate(metnet.rxns)
        lb = metnet.lb[i]
        ub = metnet.ub[i]
        lb > ub && (printstyled("rxn($i): ($rxn), lb ($lb) > ub ($ub)", 
            "\n", color = error_color); line_count += 1)
        lb > 0.0 && (printstyled("rxn($i): ($rxn), lb ($lb) > 0.0", 
            "\n", color = warn_color); line_count += 1)
        ub < 0.0 && (printstyled("rxn($i): ($rxn), ub ($ub) < 0.0", 
            "\n", color = warn_color);  line_count += 1)
        lb == ub && (printstyled("rxn($i): ($rxn), lb ($lb) == ub ($ub)", 
            "\n", color = info_color); line_count += 1)
        # TODO implement correctly the rev array
        # (isrev(metnet, i) ⊻ metnet.rev[i]) && 
        #         (printstyled("rxn($i): ($rxn), rev and bounds missmatch", "\n", 
        #                 color = error_color); line_count += 1)

        if line_count > PRINT_MAX
            printstyled("PRINT_MAX $PRINT_MAX reached!!! ... ", "\n", color = warn_color)
            break;
        end
        flush(stdout)
    end


    revscount(metnet) > 0 && printstyled("revscount: $(revscount(metnet))", "\n", color = error_color)
    fwdscount(metnet) > 0 && printstyled("fwds: $(fwdscount(metnet))", "\n", color = info_color)
    bkwdscount(metnet) > 0 && printstyled("bkwds: $(bkwdscount(metnet))", "\n", color = warn_color)
    blockscount(metnet) > 0 && printstyled("blocks: $(blockscount(metnet))", "\n", color = warn_color)
    fixxedscount(metnet) > 0 && printstyled("fixxed: $(fixxedscount(metnet))", "\n", color = info_color)

    return nothing
end