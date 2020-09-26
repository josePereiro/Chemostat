function fva(S, b, lb, ub, idxs = eachindex(lb); 
        check_obj = nothing, check_obj_atol = 1e-4,
        verbose = true, 
        upfrec = 10, 
        zeroth = 1e-10,
        on_empty_sol = (idx, sense) -> error("FBA failed, empty solution returned!!!"))

    fvalb, fvaub = (lb, ub) .|> copy
    
    sv = zeros(size(S, 2));

    if !isnothing(check_obj)
        obj_val = fba(S, b, lb, ub, check_obj).obj_val
    end
    
    n = length(idxs)
    prog = Progress(n; desc = "Doing FVA  ")
    for (i, idx) in enumerate(idxs)
        
        verbose && update!(prog, i)
        for (fva_col, sense) in [(fvalb, 1), (fvaub, -1)]
            
            sv[idx] = sense
            sol = linprog(
                sv, # Opt sense vector 
                S, # Stoichiometric matrix
                b, # row lb
                b, # row ub
                fvalb, # column lb
                fvaub, # column ub
                ClpSolver());
            x = isempty(sol.sol) ? on_empty_sol(idx, sense) : sol.sol[idx]
            fva_col[idx] = abs(x) < zeroth ? zero(x) : x
            sv[idx] = zero(sense)
        end

        # check obj
        # This just check that the obj_val is unchanged
        if !isnothing(check_obj)
            # check both first (this use the premize that only a few rxns will affect the biomass)
            fbaout = fba(S, b, fvalb, fvaub, check_obj);
            if !isapprox(fbaout.obj_val, obj_val; atol = check_obj_atol)
                # Check lb effect
                ub_ = fvaub[idx] # temp ub
                fvaub[idx] = ub[idx] # back to original ub
                fbaout = fba(S, b, fvalb, fvaub, check_obj);
                if !isapprox(fbaout.obj_val, obj_val; atol = check_obj_atol)
                    fvalb[idx] = lb[idx] # if change, back to origin
                end
                fvaub[idx] = ub_ # restore ub

                # Check ub effect
                lb_ = fvalb[idx] # temp lb
                fvalb[idx] = lb[idx] # back to original lb
                fbaout = fba(S, b, fvalb, fvaub, check_obj);
                if !isapprox(fbaout.obj_val, obj_val; atol = check_obj_atol)
                    fvaub[idx] = ub[idx] # if change, back to origin
                end
                fvalb[idx] = lb_ # restore lb
            end
        end
    end

    verbose && finish!(prog)

    return fvalb[idxs], fvaub[idxs]
end

function fva(model::MetNet, iders = eachindex(model.lb); 
        check_obj = nothing, kwargs...) 
    obj_idx = isnothing(check_obj) ? nothing : rxnindex(model, check_obj)
    idxs = [rxnindex(model, idx) for idx in iders]
    return fva(model.S, model.b, model.lb, model.ub, idxs; 
        check_obj = obj_idx, kwargs...);
end

fva(model::MetNet, ider::AbstractString; kwargs...) = 
    fva(model, [ider]; kwargs...)