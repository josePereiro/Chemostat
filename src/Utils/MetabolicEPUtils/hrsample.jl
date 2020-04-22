hrsample(metnet::MetNet; drop_samples = true, kwargs...) = 
    HRout(hrsample(metnet.S, metnet.b, metnet.lb, metnet.ub; kwargs...), drop_samples)