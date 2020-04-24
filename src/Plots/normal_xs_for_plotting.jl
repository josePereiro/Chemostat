"""
    return a number of xs points that are representative of the normal
    distribution
"""
function normal_xs_for_plotting(μ, σ, lb, ub; N = 1000)

    dx = logspace(min(log(abs(σ)), -3), 
        log(abs(ub - lb)) + 1, N)
    left_xs = μ .- dx
    right_xs = μ .+ dx
    lin_xs = collect(lb:abs(ub - lb)/N:ub)
    un = sort(unique([left_xs; right_xs; lin_xs]))
    xs = filter(x -> lb <= x <= ub, un)
    return xs
end