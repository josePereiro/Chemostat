function to_symbol_dict(src_dict)
    dict = Dict()
    for (k, dat) in src_dict
        dict[Symbol(k)] = dat
    end
    return dict
end

function struct_to_dict(obj::T) where T
    dict = Dict()
    for f in fieldnames(T)
        dict[f] = getfield(obj, f)
    end
    return dict
end