# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

struct EPFields{T<:AbstractFloat}
    av::Vector{T}
    va::Vector{T}
    a::Vector{T}
    b::Vector{T}
    μ::Vector{T}
    s::Vector{T}
    siteflagave::BitArray{1}
    siteflagvar::BitArray{1}
end

function EPFields(N::Int, expval, scalefact)
    
    siteflagvar = trues(N)
    siteflagave = trues(N)
    
    expave, expvar = parseexpval!(expval, siteflagave, siteflagvar, scalefact)
    T = typeof(scalefact)
    av = zeros(T,N)
    var = zeros(T,N)

    for (k,v) in expave
        av[k] = v
    end    
    for (k,v) in expvar
        var[k] = v
    end

    EPFields(av, var, zeros(T,N), ones(T,N), zeros(T,N), ones(T,N), siteflagave, siteflagvar)
end

EPFields(N::Int,expval::Nothing,scalefact,T) = 
    EPFields(zeros(T,N), zeros(T,N), zeros(T,N), ones(T,N), zeros(T,N), ones(T,N), trues(N), trues(N))