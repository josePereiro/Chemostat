# TODO implement correctly the rev array
update_rev!(metnet::MetNet) = nothing #metnet.rev .= ((metnet.lb .< 0.0) .& (metnet.ub .> 0.0))