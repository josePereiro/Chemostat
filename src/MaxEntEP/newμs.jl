# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function newμs(Σ,a,b,v,lb,ub,minvar,maxvar)

    Σ == 0 && (Σ = minvar)
    #lΣ = clamp(Σ,minvar,maxvar)
    #s = Σ > 0 ? clamp(inv(1.0/Σ - 1.0/b),minvar,maxvar) : minvar
    s = clamp(inv(1.0/Σ - 1.0/b), minvar,maxvar)
    μ = if Σ != b
        s * (v/Σ - a/b)
    else
        #@warn("I'm here: ub = ",ub," lb = ",lb, " Σ = ", Σ)
        0.5 * (ub+lb)
    end
    return μ,s
end