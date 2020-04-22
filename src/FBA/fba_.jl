function fba(model::MetNet, obj_ider; sense = -1.0)
    obj_idx = rxnindex(model, obj_ider)
    ei = zeros(size(model, 2));
    ei[obj_idx] = sense
    sol = linprog(
        ei, # Opt sense vector 
        model.S, # Stoichiometric matrix
        model.b, # row lb
        model.b, # row ub
        model.lb, # column lb
        model.ub, # column ub
        ClpSolver());
    return FBAout(sol.sol, sol.sol[obj_idx], obj_ider, sol)
end