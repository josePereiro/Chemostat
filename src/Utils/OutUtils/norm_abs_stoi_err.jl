norm_abs_stoi_err(metnet::MetNet, out::AbstractOut) = abs.(stoi_err(metnet.S, av(out), metnet.b)) ./ mean(abs.(av(out)))

norm_abs_stoi_err(metnet::MetNet, out::AbstractOut, ider::IDER_TYPE) = 
    norm_abs_stoi_err(metnet, out)[metindex(metnet, ider)]

