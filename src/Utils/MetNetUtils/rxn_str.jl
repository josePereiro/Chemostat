function rxn_str(metnet::MetNet, ider::IDER_TYPE)
    ridx = rxnindex(metnet, ider)
    arrow_str = isblock(metnet, ridx) ? " >< " : 
                isbkwd_bounded(metnet, ridx) ? " <== " :
                isfwd_bounded(metnet, ridx) ? " ==> " : " <==> " 
    
    reacts = rxn_reacts(metnet, ridx)
    react_str = join([string("(", S(metnet, react, ridx), ") ", 
        mets(metnet, react))  for react in reacts], " + ")
    
    prods = rxn_prods(metnet, ridx)
    prods_str = join([string("(", S(metnet, prod, ridx), ") ", 
        mets(metnet, prod))  for prod in prods], " + ")
    return react_str * arrow_str * prods_str
end
rxn_str(metnet::MetNet, iders::Vector) = [rxn_str(metnet, ider) for ider in iders]
rxn_str(metnet::MetNet) = rxn_str(metnet, metnet.rxns)