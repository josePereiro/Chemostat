function preprocess(metnet::MetNet)
    S_, b_, lb_, ub_, rxns_ = preprocess(metnet.S, metnet.b, metnet.lb, metnet.ub, metnet.rxns);
    M_, N_ = size(S_)
    MetNet(N_, M_, S_, metnet.b, metnet.c, lb_, ub_, 
        metnet.genes, metnet.rxnGeneMat, metnet.grRules, metnet.mets, rxns_, 
        metnet.metNames, metnet.metFormulas, metnet.rxnNames, metnet.rev, metnet.subSystems)
    end