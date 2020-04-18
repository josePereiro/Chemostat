hrsample(metnet::MetNet; kwargs...) = hrsample(metnet.S, metnet.b, 
                        metnet.lb, metnet.ub; kwargs...)