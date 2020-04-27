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
