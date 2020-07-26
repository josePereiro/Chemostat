function fva_check(S, b, lb, ub, idxs = eachindex(lb);
    check_obj = nothing, check_obj_atol = 1e-4)

lb_, ub_ = (lb, ub) .|> copy

sv = zeros(size(S, 2));

if !isnothing(check_obj)
    obj_val = fba(S, b, lb, ub, check_obj).obj_val
end

for (col_, col, sense) in [(lb_, lb, 1), (ub_, ub, -1)]
    
    for (i, idx) in idxs |> enumerate
        
        # minimize i
        sv[idx] = sense
        sol = linprog(
            sv, # Opt sense vector 
            S, # Stoichiometric matrix
            b, # row lb
            b, # row ub
            lb_, # column lb
            ub_, # column ub
            ClpSolver());
        isempty(sol.sol) && error("FBA failed, empty solution returned!!!")

        col_[idx] = sol.sol[idx]
        if !isnothing(check_obj)
            fbaout = fba(S, b, lb_, ub_, check_obj);
            if !isapprox(fbaout.obj_val, obj_val; atol = check_obj_atol)
                col_[idx] = col[idx]
            end
        end
    end
end

return (lb_[idxs], ub_[idxs])
end

function fva_check(model, iders = eachindex(model.rxns); check_obj = nothing, kwargs...)
obj_idx = isnothing(check_obj) ? nothing : rxnindex(model, check_obj)
idxs = [rxnindex(model, ider) for ider in iders]
fva_check(model.S, model.b, model.lb, model.ub, iders; check_obj = obj_idx, kwargs...)
end