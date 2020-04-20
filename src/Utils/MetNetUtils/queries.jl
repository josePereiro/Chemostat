isrev(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
            metnet.lb[indx] < 0.0 && metnet.ub[indx] > 0.0)

isblock(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] == 0.0 && metnet.ub[indx] == 0.0)

isfwd(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] >= 0.0 && metnet.ub[indx] > 0.0)

isbkwd(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] < 0.0 && metnet.ub[indx] <= 0.0)

isfixxed(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
    metnet.lb[indx] == metnet.ub[indx] != 0.0)

revs(metnet::MetNet) = findall((metnet.lb .< 0.0) .& (metnet.ub .> 0.0))
revscount(metnet::MetNet) = length(revs(metnet))

blocks(metnet::MetNet) = findall((metnet.lb .== 0.0) .& (metnet.ub .== 0.0))
blockscount(metnet::MetNet) = length(blocks(metnet))

fwds(metnet::MetNet) = findall((metnet.lb .>= 0.0) .& (metnet.ub .> 0.0))
fwdscount(metnet::MetNet) = length(fwds(metnet))

bkwds(metnet::MetNet) = findall((metnet.lb .< 0.0) .& (metnet.ub .<= 0.0))
bkwdscount(metnet::MetNet) = length(bkwds(metnet))

fixxeds(metnet::MetNet) = findall((metnet.lb .== metnet.ub .!= 0.0))
fixxedscount(metnet::MetNet) = length(fixxeds(metnet))

allfwd(metnet::MetNet) = fwdscount(metnet) == rxnscount(metnet)