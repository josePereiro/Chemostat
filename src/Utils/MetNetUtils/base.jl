import Base: size
size(metnet::MetNet) = size(metnet.S)
size(metnet::MetNet, dim) = size(metnet.S, dim)

