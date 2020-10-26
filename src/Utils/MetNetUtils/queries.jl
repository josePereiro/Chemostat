isrev(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
            metnet.lb[indx] < 0.0 && metnet.ub[indx] > 0.0)

isblock(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] == 0.0 && metnet.ub[indx] == 0.0)

isopen(metnet::MetNet, ider::IDER_TYPE) = !isblock(metnet::MetNet, ider::IDER_TYPE)

isfwd_bounded(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] >= 0.0 && metnet.ub[indx] > 0.0)

isbkwd_bounded(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] < 0.0 && metnet.ub[indx] <= 0.0)

isfwd_defined(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
    length(rxn_reacts(metnet, indx)) > 0)

isbkwd_defined(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
    length(rxn_prods(metnet, indx)) > 0) 

isfixxed(metnet::MetNet, ider::IDER_TYPE) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] == metnet.ub[indx] != 0.0)

revs(metnet::MetNet) = findall((metnet.lb .< 0.0) .& (metnet.ub .> 0.0))
revscount(metnet::MetNet) = length(revs(metnet))

blocks(metnet::MetNet) = findall((metnet.lb .== 0.0) .& (metnet.ub .== 0.0))
blockscount(metnet::MetNet) = length(blocks(metnet))

fwds_bounded(metnet::MetNet) = findall((metnet.lb .>= 0.0) .& (metnet.ub .> 0.0))
fwds_boundedcount(metnet::MetNet) = length(fwds_bounded(metnet))

bkwds_bounded(metnet::MetNet) = findall((metnet.lb .< 0.0) .& (metnet.ub .<= 0.0))
bkwds_boundedcount(metnet::MetNet) = length(bkwds_bounded(metnet))

fixxeds(metnet::MetNet) = findall((metnet.lb .== metnet.ub .!= 0.0))
fixxedscount(metnet::MetNet) = length(fixxeds(metnet))

allfwd(metnet::MetNet) = fwds_boundedcount(metnet) == rxnscount(metnet)

function is_exchange(metnet::MetNet, ider::IDER_TYPE)
    reacts = rxn_reacts(metnet, ider)
    prods = rxn_prods(metnet, ider)
    return xor(isempty(reacts), isempty(prods))
end
exchanges(metnet::MetNet) = findall(x -> is_exchange(metnet, x), metnet.rxns)
exchangescount(metnet::MetNet) = length(exchanges(metnet))