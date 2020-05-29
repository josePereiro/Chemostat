stoi_err(S, v, b) = S * v - b

stoi_err(metnet::MetNet, out) = stoi_err(metnet.S, av(out), metnet.b)
stoi_err(metnet::MetNet, out, ider) = 
    stoi_err(metnet.S, av(out), metnet.b)[metindex(metnet, ider)]

