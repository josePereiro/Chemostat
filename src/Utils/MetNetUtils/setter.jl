ub!(metnet::MetNet, rxnider, ub) = metnet.ub[rxnindex(metnet, rxnider)] = Float64(ub)
export ub!

lb!(metnet::MetNet, rxnider, lb) = metnet.lb[rxnindex(metnet, rxnider)] = Float64(lb)
export lb!

S!(metnet::MetNet, metider, rxnider, s) = 
    metnet.S[metindex(metnet, metider), rxnindex(metnet, rxnider)] = s

b!(metnet::MetNet, metider, b) = metnet.b[metindex(metnet, metider)] = b