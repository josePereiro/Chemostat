const warn_color = :yellow
const info_color = :blue
const error_color = :red

function _print_summary_head()
    print("SUMMARY (color code: ")
    printstyled("warning", color = warn_color)
    printstyled(", info", color = info_color)
    printstyled(", error", color = error_color)
    println(")")
end

"""
    Print useful info about the given metnet.
"""
function summary(metnet::MetNet)
    _print_summary_head()
    printstyled("model size: $(size(metnet))", "\n", color = info_color)
    _summary_bound_state(metnet)
end

function _summary_bound_state(metnet::MetNet)
    PRINT_MAX = 50 # TODO make this global?
    M, N = size(metnet)

    all(iszero(metnet.lb)) && all(iszero(metnet.ub)) && 
        (printstyled("lb and ub boths has only zero elements", "\n", color = error_color); return)
    
    # single checks
    for (name, col) in zip(["lb", "ub"], [metnet.lb, metnet.ub])
        _print_col_summary(col, name; expected_l = N, PRINT_MAX = PRINT_MAX)
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
            "\n", color = warn_color); line_count += 1)
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


    revscount(metnet) > 0 && printstyled("revscount: $(revscount(metnet))", "\n", color = warn_color)
    fwdscount(metnet) > 0 && printstyled("fwds: $(fwdscount(metnet))", "\n", color = info_color)
    bkwdscount(metnet) > 0 && printstyled("bkwds: $(bkwdscount(metnet))", "\n", color = warn_color)
    blockscount(metnet) > 0 && printstyled("blocks: $(blockscount(metnet))", "\n", color = warn_color)
    fixxedscount(metnet) > 0 && printstyled("fixxed: $(fixxedscount(metnet))", "\n", color = info_color)

    return nothing
end

function _print_col_summary(col, name; 
        expected_l = length(col), 
        PRINT_MAX = 50)

        length(col) != expected_l && (printstyled(
            " $name: ($(length(col))) != N ($expected_l), dimention missmatch", 
            "\n", color = error_color); return)
    
        unique_ = sort!(unique(col))
        length(unique_) < PRINT_MAX ?
            printstyled(" $name: $(length(unique_)) unique elment(s): ", unique_, "\n", color = info_color) :
            printstyled(" $name: $(length(unique_)) unique elment(s): min: ", 
                first(unique_), " mean: ", mean(col),
                " max: ", last(unique_), "\n", color = info_color)

end

function _print_rxn_summary(metnet, ider)
    idx = rxnindex(metnet, ider)
    printstyled(" rxn[$idx]: ", metnet.rxns[idx], " (", metnet.rxnNames[idx], ")\n", color = info_color)
    printstyled(" lb: ", metnet.lb[idx], ", ub: ", metnet.ub[idx], "\n" , color = info_color)
    printstyled(" ", rxn_str(metnet, idx), "\n" , color = info_color)
end

function _print_met_summary(metnet, ider)
    idx = metindex(metnet, ider)
    printstyled(" met[$idx]: ", metnet.mets[idx], " (", metnet.metNames[idx], ")\n", color = info_color)
    printstyled(" ", balance_str(metnet, ider), color = info_color)
end

function summary(model::MetNet, ider::IDER_TYPE)
    _print_summary_head()
    try
        _print_rxn_summary(model, ider)
    catch end
    try
        _print_met_summary(model, ider)
    catch end
end

function summary(model::MetNet, out::AbstractOut;
        PRINT_MAX = 50, 
        digits = 4)
    _print_summary_head()
    _print_out_head(out)
    _print_out_stats(model, out, PRINT_MAX)
    _print_out_exchanges_info(model, out, digits)
end

function _print_out_stats(model, out, PRINT_MAX)
    println("Out stats:")
    av_ = av(out)
    va_ = va(out)
    for (name, col) in zip(["av", "va"], [av_, va_])
        _print_col_summary(col, name; expected_l = size(model, 2), 
            PRINT_MAX = PRINT_MAX)
    end
end

function _print_out_head(fbaout::FBAout)
    printstyled(" FBA: objective val: $(fbaout.obj_val)\n", color = info_color)
end

function _print_out_head(out)
    println(" Not Implemented!!!")
end

function _print_out_exchanges_info(model::MetNet, out::AbstractOut, digits = 3)
    println("Exchange info:")
    exchs_ = model.rxns[exchanges(model)] |> sort
    for sense in [1, -1]
        for ider in exchs_
            av_ = av(model, out, ider)
            if sense * av_ > 0
                av_ = round(av_, digits = digits)
                va_ = round(va(model, out, ider), digits = digits)
                rstr = rxn_str(model, ider)
                rbounds = round.(bounds(model, ider), digits = digits)
                printstyled(" $ider => av: $av_, va: $va_, bounds: $(rbounds), eq: $(rstr)\n", color = info_color)
            end
        end
        println()
    end
end