function reaction_str(metnet, ider::Union{AbstractString, Int})
    ridx = rxnindex(metnet, ider)
    arrow_str = isblock(metnet, ridx) ? " >< " : 
                isbkwd(metnet, ridx) ? " <== " :
                isfwd(metnet, ridx) ? " ==> " : " <==> " 
    
    reacts = rxn_reacts(metnet, ridx)
    react_str = join([string("(", S(metnet, react, ridx), ") ", 
        mets(metnet, react))  for react in reacts], " + ")
    
    prods = rxn_prods(metnet, ridx)
    prods_str = join([string("(", S(metnet, prod, ridx), ") ", 
        mets(metnet, prod))  for prod in prods], " + ")
    return react_str * arrow_str * prods_str
end
reaction_str(metnet, iders) = [reaction_str(metnet, ider) for ider in iders]
reaction_str(metnet) = reaction_str(metnet, metnet.rxns)