sparsity(col::AbstractArray{<:Number}) = 1.0 - count(!iszero, col) / length(col)

function compress_dict(dict::Dict; sparsity_th = 0.66)
    new_dict = Dict()
    for (k, dat) in dict
        if dat isa AbstractArray{<:Number}
            new_dict[k] = sparsity(dat) > sparsity_th ? sparse(dat) : copy(dat)
        else
            new_dict[k] = deepcopy(dat)
        end
    end
    return new_dict
end

function uncompress_dict(dict::Dict)
    new_dict = Dict()
    for (k, dat) in dict
        if dat isa AbstractArray{<:Number} && issparse(dat)
            new_dict[k] = dat |> collect
        else
            new_dict[k] = deepcopy(dat)
        end
    end
    return new_dict
end