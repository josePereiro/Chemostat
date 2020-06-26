# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
function fva_preprocess(S,b,lb,ub,rxns; 
        verbose = false, eps = 0.0, 
        ignored = [], # rxns skip preprocess
        upfrec = 10)
        
    b, lb, ub = copy(b), copy(lb), copy(ub)
    n = length(lb)
    ei = zeros(n)
    blocked = falses(n)
    upfrec = floor(Int, n/upfrec)
    for i=1:n

        i in ignored && continue
        
        show_progress = verbose && (i == 1 || i % upfrec == 0 || i == n)
        show_progress && (print("fva_processing [$i / $n]        \r"); flush(stdout))

        ei[i] = -1.0
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
    verbose && (println("fva_processing done!!!                                           "); flush(stdout))
    idxb = findall(blocked)

    verbose && println("$(sum(blocked)) blocked fluxes")

    unblocked = (lb .< ub)

    return S[:,unblocked], b - S*(blocked .* lb), lb[unblocked], ub[unblocked], rxns[unblocked], idxb
end


function fva_preprocess(metnet::MetNet; 
            verbose = false, eps = 0.0, return_blocked = false,
            ignored = [] # rxns skip preprocess
        )
    ignored = map((r) -> rxnindex(metnet, r), ignored)

    S_, b_, lb_, ub_, rxns_, idxb_ = fva_preprocess(metnet.S, metnet.b, metnet.lb, 
        metnet.ub, metnet.rxns; verbose = verbose, eps = eps, ignored = ignored);
    
    metnet = MetNet(metnet; S = S_, b = b_, lb = lb_, ub = ub_, rxns = rxns_)
    return return_blocked ? (metnet, idxb_) : metnet
end