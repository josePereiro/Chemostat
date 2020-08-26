#TODO make test for this
"""
    Example of 'intake_info' = Dict("gt" => Dict("ub" => 100.0, "c" => 10.0))
"""
function apply_bound!(metnet::MetNet, ξ, intake_info::Dict = Dict();
        emptyfirst = false, ignore_miss = false)
    ξ < 0 && error("ξ must be positive")

    emptyfirst && empty!(metnet.intake_info)
    for (exch, exch_info) in intake_info
        ignore_miss && !(exch in metnet.rxns) && continue
        metnet.intake_info[exch] = exch_info
    end
    
    for (exch, exch_info) in metnet.intake_info
        
        !haskey(exch_info, "c") && error("'c' not defined for exchange '$(rxns(metnet, exch))'");
        
        c = exch_info["c"]
        c < 0.0 && error("negative 'c' ($c) in exchange '$(rxns(metnet, exch))'")
        
        if haskey(exch_info, "ub")
            ub = exch_info["ub"]
            ub!(metnet, exch, stst_bound(c, ξ, ub))
        elseif haskey(exch_info, "lb")
            lb = exch_info["lb"]
            lb!(metnet, exch, stst_bound(c, ξ, lb))
        else
            error("neither 'lb' or 'ub' defined for exchange '$(rxns(metnet, exch))'")
        end
    end
    return metnet
end