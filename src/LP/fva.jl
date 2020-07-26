function fva(S, b, lb, ub, idxs = eachindex(lb); 
        check_obj = nothing, check_obj_atol = 1e-4,
        verbose = true, 
        upfrec = 10, 
        zeroth = 1e-10)

    fvalb, fvaub = (lb, ub) .|> copy
    
    sv = zeros(size(S, 2));

    if !isnothing(check_obj)
        obj_val = fba(S, b, lb, ub, check_obj).obj_val
    end
    
    n = length(idxs)
    for (fva_col, orig_col, sense) in [(fvalb, lb, 1), (fvaub, ub, -1)]

        for (i, idx) in enumerate(idxs)

            show_progress = verbose && sense == 1 && (i == 1 || i % upfrec == 0 || i == n)
            show_progress && (print("fva[$i / $n]        \r"); flush(stdout))

            sv[idx] = sense
            sol = linprog(
                sv, # Opt sense vector 
                S, # Stoichiometric matrix
                b, # row lb
                b, # row ub
                fvalb, # column lb
                fvaub, # column ub
                ClpSolver());
            isempty(sol.sol) && error("FBA failed, empty solution returned!!!")
            
            fva_col[idx] = sol.sol[idx]
            if !isnothing(check_obj)
                fbaout = fba(S, b, fvalb, fvaub, check_obj);
                if !isapprox(fbaout.obj_val, obj_val; atol = check_obj_atol)
                    fva_col[idx] = orig_col[idx]
                end
            end
            sv[idx] = zero(sense)
        end
    end

    verbose && (println("done!!!", " "^100); flush(stdout))
    
    return fvalb[idxs], fvaub[idxs]
end

function fva(model::MetNet, iders = eachindex(model.lb); 
        check_obj = nothing, check_obj_atol = 1e-4,
        verbose = true, upfrec = 10) 
    obj_idx = isnothing(check_obj) ? nothing : rxnindex(model, check_obj)
    idxs = [rxnindex(model, idx) for idx in iders]
    return fva(model.S, model.b, model.lb, model.ub, idxs; 
        check_obj = obj_idx, check_obj_atol = check_obj_atol,
        verbose = verbose, upfrec = upfrec);
end

fva(model::MetNet, ider::AbstractString; kwargs...) = 
    fva(model, [ider]; kwargs...)