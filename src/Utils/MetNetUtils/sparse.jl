function SparseArrays.sparse(model::MetNet, fields::Vector{Symbol} = [:S])
    net = Dict()
    for f in fields
        net[f] = getfield(model, f) |> sparse
    end
    return MetNet(model; reshape = false, net...)
end
