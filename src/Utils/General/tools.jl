function to_symbol_dict(src_dict)
    dict = Dict()
    for (k, dat) in src_dict
        dict[Symbol(k)] = dat
    end
    return dict
end