# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function eponesweep!(epfields::EPFields{T}, epalg::EPAlg, epmat::EPMat, stat) where T
    @extract epfields : av va a b μ s siteflagave siteflagvar
    @extract epalg : alpha beta_vec minvar maxvar epsconv damp
    @extract epmat : KKc dKKbk invKKPD lb ub KY v



    # KK = βF'F
    # KK = KKc if diag(KKc) = dKKbk
    # Σ^-1 = (KK + D)
    KKdi = diagind(KKc);
    KKc[KKdi] = dKKbk + (1.0 ./ b)
    
    # Σ
    stat[:elapsed_eponesweep_inv] = @elapsed inplaceinverse!(invKKPD,KKc)

    
    
    # KY = βF'y
    # v¯ = Σ(KY + Da) (original ep)
    # mul!(v,invKKPD, (KY + a./b))
    # v¯ = Σ(KY + Da) + Σ * beta_vec (maxent-ep)
    v .= invKKPD * (KY + a./b) + invKKPD * beta_vec
    

    errav = errva = errμ = errs = typemin(T)
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