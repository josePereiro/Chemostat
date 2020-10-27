function clampfields!(obj, fields::Vector{Symbol}; abs_max::Real = 1e3, zeroth::Real = 1e-8)
    for f in fields
        dat = getfield(obj, f)
        @inbounds for i in eachindex(dat)
            datum = dat[i]
            dat[i] = clamp(abs(datum) < zeroth ? zero(datum) : datum, -abs_max, abs_max)
        end
    end
    return obj
end