isrev(metnet::MetNet, ider) = (indx = rxnindex(metnet, ider); 
            metnet.lb[indx] < 0.0 && metnet.ub[indx] > 0.0)