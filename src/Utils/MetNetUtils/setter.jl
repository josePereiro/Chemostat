ub!(metnet::MetNet, rxnider, ub::Real) = metnet.ub[rxnindex(metnet, rxnider)] = Float64(ub)
export ub!

lb!(metnet::MetNet, rxnider, lb::Real) = metnet.lb[rxnindex(metnet, rxnider)] = Float64(lb)
export lb!

bounds!(metnet::MetNet, ider, lb::Real, ub::Real) = 
    (idx = rxnindex(metnet, ider); metnet.lb[idx] = lb; metnet.ub[idx] = ub; nothing)

S!(metnet::MetNet, metider, rxnider, s::Real) = 
    metnet.S[metindex(metnet, metider), rxnindex(metnet, rxnider)] = s

b!(metnet::MetNet, metider, b) = metnet.b[metindex(metnet, metider)] = b