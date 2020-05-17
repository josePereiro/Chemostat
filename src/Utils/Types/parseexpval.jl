# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function parseexpval!(expval, siteflagave::BitArray{1}, siteflagvar::BitArray{1},scalefact::Float64)

    expave, expvar=_parseexpval!(expval,siteflagave,siteflagvar)

    for (k,v) in expave
        expave[k] = v/scalefact
    end

    for (k,v) in expvar
        expvar[k] = v/scalefact^2
    end
    expave, expvar
end

function _parseexpval!(expval::Tuple,siteflagave::BitArray{1},siteflagvar::BitArray{1})

    N = length(siteflagave)
    length(expval) == 3 || error("We expect a 3-uple here")
    expsite = expval[3]
    1 <= expsite <= N || error("expsite = $expsite not ∈ 1,...,$N")
    expave = Dict{Int,Float64}()
    expvar = Dict{Int,Float64}()

    if expval[1] != nothing
        siteflagave[expsite] = false
        expave[expsite] = expval[1]
    end
    if expval[2] != nothing
        siteflagvar[expsite] = false
        expvar[expsite] = expval[2]
    end
    expave,expvar
end

function _parseexpval!(expval::Vector,siteflagave::BitArray{1},siteflagvar::BitArray{1})

    N = length(siteflagave)
    expave = Dict{Int,Float64}()
    expvar = Dict{Int,Float64}()
    for i in eachindex(expval)
        expsite = expval[i][3]
        1 <= expsite <= N || error("expsite = $minsite not ∈ 1,...,$N")
        if expval[i][1] != nothing
            siteflagave[expsite] = false
            expave[expsite] = expval[i][1]
        end
        if expval[i][2] != nothing
            siteflagvar[expsite] = false
            expvar[expsite] =  expval[i][2]
        end
    end

    expave,expvar
end
_parseexpval!(nothing,siteflagave::BitArray{1},siteflagvar::BitArray{1})=(Dict{Int,Float64}(),Dict{Int,Float64}())