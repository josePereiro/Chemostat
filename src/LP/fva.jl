function fva(S, b, lb, ub, idxs = eachindex(lb); 
        verbose = true, 
        upfrec = 10, 
        zeroth = 1e-10)

    lvals = zeros(eltype(S), length(idxs))
    uvals = zeros(eltype(S), length(idxs))
    
    sv = zeros(size(S, 2));
    
    n = length(idxs)
    for (col, sense) in [(lvals, 1), (uvals, -1)]

        for (i, idx) in enumerate(idxs)

            show_progress = verbose && sense == 1 && (i == 1 || i % upfrec == 0 || i == n)
            show_progress && (print("fva[$i / $n]        \r"); flush(stdout))

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
            
            x = sol.sol[idx]
            col[i] = abs(x) < zeroth ? zero(x) : x
            sv[idx] = zero(sense)
        end
    end

    verbose && (println("done!!!", " "^100); flush(stdout))
    
    return lvals, uvals
end

function fva(model::MetNet, iders = eachindex(model.lb); verbose = true, upfrec = 10) 
    idxs = [rxnindex(model, idx) for idx in iders]
    return fva(model.S, model.b, model.lb, model.ub, idxs; 
        verbose = verbose, upfrec = upfrec);
end

fva(model, ider::AbstractString; verbose = true, upfrec = 10) = 
    fva(model, [ider]; verbose = verbose, upfrec = upfrec)