S!(metnet::MetNet, metider::IDER_TYPE, rxnider::IDER_TYPE, s::Real) =   
    metnet.S[metindex(metnet, metider), rxnindex(metnet, rxnider)] = s

b!(metnet::MetNet, metider::IDER_TYPE, b::Real) = metnet.b[metindex(metnet, metider)] = b

ub!(metnet::MetNet, rxnider::IDER_TYPE, ub::Real) = metnet.ub[rxnindex(metnet, rxnider)] = Float64(ub)

lb!(metnet::MetNet, rxnider::IDER_TYPE, lb::Real) = metnet.lb[rxnindex(metnet, rxnider)] = Float64(lb)

bounds!(metnet::MetNet, ider::IDER_TYPE, lb::Real, ub::Real) = 
    (idx = rxnindex(metnet, ider); metnet.lb[idx] = lb; metnet.ub[idx] = ub; nothing)

