isrev(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
            metnet.lb[indx] < 0.0 && metnet.ub[indx] > 0.0)

revs(metnet::MetNet) = findall((metnet.lb .< 0.0) .& (metnet.ub .> 0.0))
revscount(metnet::MetNet) = length(revs(metnet))