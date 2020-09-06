# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function eponesweepT0!(epfields::EPFields, epalg::EPAlg, epmatT0::EPMatT0, stat = Dict())
    @extract epfields : av va a b μ s siteflagave siteflagvar
    @extract epalg : alpha beta_vec minvar maxvar epsconv damp
    @extract epmatT0 : Σy Σw G lb ub vy vw Y


    M = size(G,1)
    N = length(av)

    idxy = 1:M # dependent variables
    idxw = M+1:N # independent variables

    ay,aw = view(a,idxy), view(a,idxw) # dep an ind prior mean (epfields)
    by,bw = view(b,idxy), view(b,idxw) # dep an ind prior variance (epfields)
    sy,sw = view(s,idxy), view(s,idxw) # dep and ind untruncated marginals variance (epfields)
    μy,μw = view(μ,idxy), view(μ,idxw) # dep and ind untruncated marginals mean (epfields)
    avy,avw = view(av,idxy), view(av,idxw) # dep and ind truncated marginals mean (epfields)
    vay,vaw = view(va,idxy), view(va,idxw) # dep and ind truncated marginals variance (epfields)

    minerr = typemin(av[1])
    errav,errva,errμ,errs = minerr,minerr,minerr,minerr
    
    # All fields in epmat are updated from the epfields of last sweep
    # (?) covariance matrix of independent variables (epmat)
    stat[:elapsed_eponesweep_inv] = @elapsed begin
        Σw = inv(Diagonal(1.0 ./ bw) + G' * Diagonal( 1.0 ./ by ) * G) 
    end
    #fast_similarity_inv!(Σw, bw,  by, G)
    mul!(Σy,G*Σw,G') # (?) covariance matrix of dependent variables (epmat)
    # Original ep
    # mul!(vw,Σw, aw ./ bw - G'*(ay ./ by)) # (?) mean vector of independent variables (epmat)
    # ep-maxent
    vw = Σw * (aw ./ bw - G'*(ay ./ by)) + Σw * beta_vec

    mul!(vy,G,vw) # (?) mean vector of dependent variables (epmat)
    for i in eachindex(vy) vy[i] = -vy[i] + Y[i] end

    for i in eachindex(μw)  # loop M+1:N
        newμw,newsw = newμs(Σw[i,i],aw[i],bw[i],vw[i],lb[i+M],ub[i+M], minvar,maxvar)
        errμ = max(errμ, abs(μw[i]-newμw))
        errs = max(errs, abs(sw[i]-newsw))
        μw[i] = newμw
        sw[i] = newsw
        # println("μw[$(i+M)] = ", μw[i]," sw[$(i+M)] = ", sw[i], " Σw = ",Σw[i,i] )


        newavw,newvaw = newav(sw[i],μw[i],avw[i],vaw[i],siteflagave[i+M],siteflagvar[i+M],
                              lb[i+M],ub[i+M],minvar,maxvar)
        errav = max(errav,abs(avw[i]-newavw))
        errva = max(errva,abs(vaw[i]-newvaw))
        avw[i] = newavw
        vaw[i] = newvaw

        newaw,newbw = matchmom(μw[i],sw[i],avw[i],vaw[i],minvar,maxvar)
        aw[i] = damp * aw[i] + (1.0-damp)*newaw # modify a in epfields
        bw[i] = damp * bw[i] + (1.0-damp)*newbw # modify b in epfields
    end

    for i in eachindex(μy)   # loop  1:M

        newμy,newsy = newμs(Σy[i,i],ay[i],by[i], vy[i],lb[i],ub[i],minvar,maxvar)
        errμ = max(errμ, abs(μy[i]-newμy))
        errs = max(errs, abs(sy[i]-newsy))
        μy[i] = newμy # modify μ in epfields
        sy[i] = newsy # modify s in epfields
#        println("μy[$i] = ", μy[i]," sy[$i] = ", sy[i], " Σ = ", Σy[i,i], " (",lb[i],":",ub[i],")"," ay[$i] = ",ay[i], " by[$i] = ", by[i])

        newavy,newvay = newav(sy[i],μy[i],avy[i],vay[i],siteflagave[i],siteflagvar[i],
                              lb[i],ub[i],minvar,maxvar)
        errav = max(errav,abs(avy[i]-newavy))
        errva = max(errva,abs(vay[i]-newvay))
        avy[i] = newavy # modify av in epfields
        vay[i] = newvay # modify va in epfields

        neway,newby = matchmom(μy[i],sy[i],avy[i],vay[i],minvar,maxvar)
        ay[i] = damp * ay[i] + (1.0-damp)*neway # modify a in epfields
        by[i] = damp * by[i] + (1.0-damp)*newby # modify b in epfields
    end
    return errav, errva, errμ, errs
end