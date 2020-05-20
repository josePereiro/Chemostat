# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
function preprocess(S,b,lb,ub,rxns; verbose = false, eps = 0.0)
    b,lb, ub = copy(b), copy(lb), copy(ub)
    n = length(lb)
    ei = zeros(n)
    blocked = falses(n)
    for i=1:n
        ei[i] = -1.0
        verbose && (print("processing [$i / $n]        \r"); flush(stdout))
        sol = linprog(ei, S, b, b, lb, ub, ClpSolver())
        isempty(sol.sol) && error("FVA fails to find a solution maximixing rxn $(rxns[i])")
        ub[i] = min(ub[i], sol.sol[i])
        ei[i] = +1.0
        sol = linprog(ei, S, b, b, lb, ub, ClpSolver())
        isempty(sol.sol) && error("FVA fails to find a solution minimaxing rxn $(rxns[i])")
        lb[i] = max(lb[i], sol.sol[i])
        ei[i] = 0.0

        if lb[i] == ub[i]
            blocked[i] = true
            ub[i] += eps
            lb[i] -= eps
        end
    end
    verbose && (println("done!!!                                           "); flush(stdout))
    idxb = findall(blocked)

    verbose && println("$(sum(blocked)) blocked fluxes")

    unblocked = (lb .< ub)

    return S[:,unblocked], b - S*(blocked .* lb), lb[unblocked], ub[unblocked], rxns[unblocked], idxb
end


function preprocess(metnet::MetNet; verbose = false, eps = 0.0, return_blocked = false,)
    S_, b_, lb_, ub_, rxns_, idxb_ = preprocess(metnet.S, metnet.b, metnet.lb, 
        metnet.ub, metnet.rxns; verbose = verbose, eps = eps);
    
    metnet = MetNet(metnet; S = S_, b = b_, lb = lb_, ub = ub_, rxns = rxns_)
    return return_blocked ? (metnet, idxb_) : metnet
end