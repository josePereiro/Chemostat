norm1_stoi_err(S, v, b) = abs.(stoi_err(S, v, b)) ./ mean(abs.(v))

norm1_stoi_err(metnet::MetNet, out::AbstractOut) = norm1_stoi_err(metnet.S, av(out), metnet.b)

norm1_stoi_err(metnet::MetNet, out::AbstractOut, ider::IDER_TYPE) = 
    norm1_stoi_err(metnet, out)[metindex(metnet, ider)]

function norm2_stoi_err(model, out; 
    xfun = av)
    errs = stoi_err(model, out)
    x = xfun(out)
    max_abs_flxs = []
    for meti in eachindex(model.mets)
        rxnis = met_rxns(model, meti)
        push!(max_abs_flxs, isempty(rxnis) ? 1.0 : maximum(abs.(x[rxnis])))
    end
    norm_errs =[]
    for (err, norm) in zip(errs, max_abs_flxs)
        norm_err = err == 0.0 ? 0.0 :
            abs(err) / norm
        push!(norm_errs, norm_err)
    end
    norm_errs
end
norm2_stoi_err(metnet::MetNet, out::AbstractOut, ider::IDER_TYPE) = 
    norm2_stoi_err(metnet, out)[metindex(metnet, ider)]