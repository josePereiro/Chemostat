# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

struct EPout{T<:AbstractFloat}
    μ::Vector{T}
    σ::Vector{T}
    av::Vector{T}
    va::Vector{T}
    sol::EPFields{T}
    status::Symbol
end