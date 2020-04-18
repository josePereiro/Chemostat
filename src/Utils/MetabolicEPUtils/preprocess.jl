preprocess(metnet::MetNet) = MetNet(
    preprocess(metnet.S, metnet.b, metnet.lb, metnet.ub, metnet.rxns)...)