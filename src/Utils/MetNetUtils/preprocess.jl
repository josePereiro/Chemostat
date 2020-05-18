# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
function preprocess(S,b,lb,ub,rxns; verbose = false)
    b,lb, ub = copy(b), copy(lb), copy(ub)
    n = length(lb)
    ei = zeros(n)
    for i=1:n
        ei[i] = -1.0
        sol = linprog(ei, S, b, b, lb, ub, ClpSolver())
        ub[i] = min(ub[i], sol.sol[i])
        ei[i] = +1.0
        sol = linprog(ei, S, b, b, lb, ub, ClpSolver())
        lb[i] = max(lb[i], sol.sol[i])
        ei[i] = 0.0
    end
    blocked = (lb .== ub)
    idxb = findall(lb .== ub)
    if verbose
        println("$(sum(blocked)) blocked fluxes")
        for i = 1:sum(blocked)
            println(rxns[idxb[i]], "is fixed to ", lb[idxb[i]])
        end
    end
    unblocked = (lb .< ub)
    return S[:,unblocked], b - S*(blocked .* lb), lb[unblocked], ub[unblocked], rxns[unblocked]
end


function preprocess(metnet::MetNet; verbose = false)
    S_, b_, lb_, ub_, rxns_ = preprocess(metnet.S, metnet.b, metnet.lb, 
        metnet.ub, metnet.rxns; verbose = verbose);
    M_, N_ = size(S_)
    return MetNet(N_, M_, S_, metnet.b, metnet.c, lb_, ub_, 
                    metnet.genes, metnet.rxnGeneMat, metnet.grRules, metnet.mets, rxns_, 
                    metnet.metNames, metnet.metFormulas, metnet.rxnNames, 
                    ((lb_ .< 0.0) .& (ub_ .> 0.0)), metnet.subSystems)
end