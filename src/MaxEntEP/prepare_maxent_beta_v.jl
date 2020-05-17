# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function prepare_maxent_beta_v(epmat::EPMat, beta_maxent::Real, obj_idx::Integer)
    N = length(epmat.v)
    beta_mv = zeros(N)
    beta_maxent == 0.0 && return beta_mv
    0 >= obj_idx > N && error("You must specified a valid 'obj_idx' ∈ [0, $N], current $obj_idx")
    beta_mv[obj_idx] = beta_maxent
    return beta_mv
end

function prepare_maxent_beta_v(epmat::EPMatT0, beta_maxent::Real, obj_idx::Integer)
        M = length(epmat.vy) # Number of dependent variables
        N = length(epmat.vy) + length(epmat.vw)
        beta_mv = zeros(N - M) # same length as independent variables
        beta_maxent == 0.0 && return beta_mv
        0 >= obj_idx > N && error("You must specified a valid 'obj_idx' ∈ [0, $N], current $obj_idx")
        obj_idx_epmat = findfirst(isequal(obj_idx), epmat.idx[M + 1:end])
        isnothing(obj_idx_epmat) && error("obj_idx must be an independen variable, try move it to the last indexes");
        beta_mv[obj_idx_epmat] = beta_maxent
end