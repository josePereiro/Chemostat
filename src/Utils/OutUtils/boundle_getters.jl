# Boundle
va_ep(boundle::ChstatBoundle, ξ, β) = 
    va(get_epout(boundle, ξ, β))
va_ep(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    va(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)

av_ep(boundle::ChstatBoundle, ξ, β) = 
    av(get_epout(boundle, ξ, β))
av_ep(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    av(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)

μ_ep(boundle::ChstatBoundle, ξ, β) = 
    μ(get_epout(boundle, ξ, β))
μ_ep(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    μ(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)

σ_ep(boundle::ChstatBoundle, ξ, β) = 
    σ(get_epout(boundle, ξ, β))
σ_ep(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    σ(get_metnet(boundle, ξ), get_epout(boundle, ξ, β), ider)

va_hr(boundle::ChstatBoundle, ξ, β) = 
    va(get_hrout(boundle, ξ, β))
va_hr(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    va(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)

av_hr(boundle::ChstatBoundle, ξ, β) = 
    av(get_hrout(boundle, ξ, β))
av_hr(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    av(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)

μ_hr(boundle::ChstatBoundle, ξ, β) = 
    μ(get_hrout(boundle, ξ, β))
μ_hr(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    μ(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)

σ_hr(boundle::ChstatBoundle, ξ, β) = 
    σ(get_hrout(boundle, ξ, β))
σ_hr(boundle::ChstatBoundle, ξ, β, ider::Union{AbstractString, Integer}) = 
    σ(get_metnet(boundle, ξ), get_hrout(boundle, ξ, β), ider)


va_fba(boundle::ChstatBoundle, ξ) = va(get_fbaout(boundle, ξ))
va_fba(boundle::ChstatBoundle, ξ, ider::Union{AbstractString, Integer}) = 
    va(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)

av_fba(boundle::ChstatBoundle, ξ) = 
    av(get_fbaout(boundle, ξ))
av_fba(boundle::ChstatBoundle, ξ, ider::Union{AbstractString, Integer}) = 
    av(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)

μ_fba(boundle::ChstatBoundle, ξ) = 
    μ(get_fbaout(boundle, ξ))
μ_fba(boundle::ChstatBoundle, ξ, ider::Union{AbstractString, Integer}) = 
    μ(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)

σ_fba(boundle::ChstatBoundle, ξ) = 
    σ(get_fbaout(boundle, ξ))
σ_fba(boundle::ChstatBoundle, ξ, ider::Union{AbstractString, Integer}) = 
    σ(get_metnet(boundle, ξ), get_fbaout(boundle, ξ), ider)