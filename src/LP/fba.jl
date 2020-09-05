function fba(S, b, lb, ub, obj_idx::Integer; sense = -1.0, 
        on_empty_sol = () -> error("FBA failed, empty solution returned!!!"))
    sv = zeros(size(S, 2));
    sv[obj_idx] = sense
    sol = linprog(
        sv, # Opt sense vector 
        S, # Stoichiometric matrix
        b, # row lb
        b, # row ub
        lb, # column lb
        ub, # column ub
        ClpSolver());
    isempty(sol.sol) && on_empty_sol()
    return FBAout(sol.sol, sol.sol[obj_idx], obj_idx, sol)
end

function fba(S, b, lb, ub, obj_idx::Integer, cost_idx::Integer;
        on_empty_sol = () -> error("FBA failed, empty solution returned!!!"))
    # maximizing obj
    fbaout1 = fba(S, b, lb, ub, obj_idx; sense = -1.0, 
        on_empty_sol = on_empty_sol)
    # fix obj
    lb_ = copy(lb)
    ub_ = copy(ub)
    ub_[obj_idx] = fbaout1.obj_val
    lb_[obj_idx] = fbaout1.obj_val
    # minimize cost
    return fba(S, b, lb_, ub_, cost_idx; sense = 1.0,
        on_empty_sol = on_empty_sol)
end

function fba(model::MetNet, obj_ider::IDER_TYPE; kwargs...)
    obj_idx = rxnindex(model, obj_ider)
    return fba(_extract_dense(model, [:S, :b, :lb, :ub])..., obj_idx; kwargs...)
end

function fba(model::MetNet, obj_ider::IDER_TYPE, cost_ider::IDER_TYPE; kwargs...)
    obj_idx = rxnindex(model, obj_ider)
    cost_idx = rxnindex(model, cost_ider)
    return fba(_extract_dense(model, [:S, :b, :lb, :ub])..., obj_idx, cost_idx; kwargs...)
end

function _extract_dense(model, fields)
    extracted = []
    for f in fields
        dat = getfield(model, f)
        if dat isa AbstractArray{<:Number} && issparse(dat)
            dat = Array(dat)
        end
        push!(extracted, dat)
    end
    return extracted
end