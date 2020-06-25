function fva(S, b, lb, ub, idxs = eachindex(lb))
    
    lvals = zeros(eltype(S), length(idxs))
    uvals = zeros(eltype(S), length(idxs))
    
    sv = zeros(size(S, 2));
    
    for (col, sense) in [(lvals, 1), (uvals, -1)]
        for (i, idx) in enumerate(idxs)
            sv[idx] = sense
            sol = linprog(
                sv, # Opt sense vector 
                S, # Stoichiometric matrix
                b, # row lb
                b, # row ub
                lb, # column lb
                ub, # column ub
                ClpSolver());
            isempty(sol.sol) && error("FBA failed, empty solution returned!!!")
            
            col[i] = sol.sol[idx]
            sv[idx] = zero(sense)
        end
    end
    
    return lvals, uvals
end

function fva(model, iders = eachindex(model.lb)) 
    idxs = [rxnindex(model, idx) for idx in iders]
    return fva(model.S, model.b, model.lb, model.ub, idxs);
end

fva(model, ider::AbstractString) = fva(model, [ider])