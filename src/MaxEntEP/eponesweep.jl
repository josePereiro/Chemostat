# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function eponesweep!(X::EPFields,epalg::EPAlg, epmat::EPMat)
    @extract X : av va a b μ s siteflagave siteflagvar
    @extract epalg : alpha beta_vec minvar maxvar epsconv damp
    @extract epmat : KK KKPD invKKPD lb ub KY v

    for i in eachindex(b)
        # KK = βF'F
        # Σ^-1 = (KK + D)
        KKPD[i,i] = KK[i,i] + 1.0/b[i]
    end
    # Σ
    inplaceinverse!(invKKPD,KKPD)

    minerr = typemin(av[1])
    errav,errva,errμ,errs = minerr,minerr,minerr,minerr
    
    # KY = βF'y
    # v¯ = Σ(KY + Da) (original ep)
    # mul!(v,invKKPD, (KY + a./b))
    # v¯ = Σ(KY + Da) + Σ * beta_vec (maxent-ep)
    v .= invKKPD * (KY + a./b) + invKKPD * beta_vec
    

    for i in eachindex(av)
        # Parameters of the tilted (Confirmed)
        newμ,news = newμs(invKKPD[i,i],a[i],b[i],v[i],lb[i],ub[i],minvar, maxvar)
        errμ = max(errμ, abs(μ[i]-newμ))
        errs = max(errs, abs(s[i]-news))
        μ[i] = newμ
        s[i] = news

        # Parameters of the truncated marginals
        newave,newva = newav(s[i],μ[i],av[i],va[i],siteflagave[i],siteflagvar[i],lb[i],ub[i],minvar,maxvar)
        errav = max(errav,abs(av[i]-newave))
        errva = max(errva,abs(va[i]-newva))
        av[i] = newave
        va[i] = newva

        # Parameter of the prior
        newa,newb = matchmom(μ[i],s[i],av[i],va[i],minvar,maxvar)
        a[i] = damp * a[i] + (1.0-damp)*newa
        b[i] = damp * b[i] + (1.0-damp)*newb
    end
    return errav,errva,errμ,errs
end