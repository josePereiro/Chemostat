# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function warmup(S, b, lb, ub)
    n = length(lb)
    x0 = zeros(n)
    ei = zeros(n)
    for i=1:n
        ei[i] = -1.0
        sol=linprog(ei, S, b, b, lb, ub, ClpSolver())
        x0 .+= sol.sol / 2n
        ei[i] = +1.0
        sol=linprog(ei, S, b, b, lb, ub, ClpSolver())
        ei[i] = 0.0
        x0 .+= sol.sol / 2n
    end
    x0
end