function check_dims(model)
    M, N = size(model)
    for (expected_len, field) in [
            (N, :lb), (N, :ub), (M, :b), 
            (N, :c), (M, :mets), (N, :rxns)
        ]
        dat = getfield(model, field)
        dat_len = length(dat)
        expected_len != dat_len && 
            error("'$(field)' len = $(dat_len) missmatch model size = $(size(model))")
    end
end
