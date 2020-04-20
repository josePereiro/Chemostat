ub!(metnet::MetNet, ider, ub) = metnet.ub[rxnindex(metnet, ider)] = Float64(ub)
export ub!

lb!(metnet::MetNet, ider, lb) = metnet.lb[rxnindex(metnet, ider)] = Float64(lb)
export lb!

S!(metnet::MetNet, metider, rxnider, s) = 
    metnet.S[metindex(metnet, metider), rxnindex(metnet, rxnider)] = s