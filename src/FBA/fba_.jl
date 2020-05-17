function fba(S, b, lb, ub, obj_idx::Integer; sense = -1.0)
    ei = zeros(size(S, 2));
    ei[obj_idx] = sense
    sol = linprog(
        ei, # Opt sense vector 
        S, # Stoichiometric matrix
        b, # row lb
        b, # row ub
        lb, # column lb
        ub, # column ub
        ClpSolver());
    return FBAout(sol.sol, sol.sol[obj_idx], obj_idx, sol)
end

function fba(model::MetNet, obj_ider; sense = -1.0)
    obj_idx = rxnindex(model, obj_ider)
    return fba(model.S, model.b, model.lb, model.ub, obj_idx; sense = sense)
end