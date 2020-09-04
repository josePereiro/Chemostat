norm_abs_stoi_err(S, v, b) = abs.(stoi_err(S, v, b)) ./ mean(abs.(v))

norm_abs_stoi_err(metnet::MetNet, out::AbstractOut) = norm_abs_stoi_err(metnet.S, av(out), metnet.b)

norm_abs_stoi_err(metnet::MetNet, out::AbstractOut, ider::IDER_TYPE) = 
    norm_abs_stoi_err(metnet, out)[metindex(metnet, ider)]

