# A struct that contain all the data required for converging ep
struct EPModel{T<:Real}
    scalefact::T
    updatealg!::Function
    epfields::EPFields{T}
    epmat::AbstractEPMat
    alpha::T
    beta_vec::SparseVector{T, Int}
end
