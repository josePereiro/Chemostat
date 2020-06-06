# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
# The original type did not include beta_maxent vector

struct EPAlg{T<:AbstractFloat}
    alpha::T
    beta_vec::Vector{T}
    minvar::T
    maxvar::T
    epsconv::T
    damp::T
    maxiter::Int
    verbose::Bool
end