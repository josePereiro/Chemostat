empty_epout(N, T = Float64) = EPout(zeros(T,N), ones(T,N), zeros(T,N), ones(T,N), EPFields(N, nothing, T), :converged, 1)

empty_epout(model::MetNet) = empty_epout(size(model, 2), eltype(model.S))