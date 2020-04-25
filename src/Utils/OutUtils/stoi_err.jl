stoi_err(S, v, b) = S * v - b

av_stoi_err(metnet::MetNet, out) = stoi_err(metnet.S, va(out), metnet.b)
av_stoi_err(metnet::MetNet, out, ider) = 
    stoi_err(metnet.S, va(out), metnet.b)[metindex(metnet, ider)]

av_stoi_err_ep(bound::ChstatBoundle, ξ::Real, β::Real) = 
    av_stoi_err(get_metnet(bound, ξ), get_epout(bound, ξ, β))
av_stoi_err_ep(bound::ChstatBoundle, ξ::Real, β::Real, ider) = 
    av_stoi_err(get_metnet(bound, ξ), get_epout(bound, ξ, β), ider)

av_stoi_err_fba(bound::ChstatBoundle, ξ::Real) = 
    av_stoi_err(get_metnet(bound, ξ), get_fbaout(bound, ξ))
av_stoi_err_fba(bound::ChstatBoundle, ξ::Real, ider) = 
    av_stoi_err(get_metnet(bound, ξ), get_fbaout(bound, ξ), ider)

av_stoi_err_hr(bound::ChstatBoundle, ξ::Real, β::Real) = 
    av_stoi_err(get_metnet(bound, ξ), get_hrout(bound, ξ, β))
av_stoi_err_hr(bound::ChstatBoundle, ξ::Real, β::Real, ider) = 
    av_stoi_err(get_metnet(bound, ξ), get_hrout(bound, ξ, β), ider)
