# A struct that contain all the data required for converging ep
struct EPModel{T<:Real}
    scalefact::T
    updatealg!::Function
    epfield::EPFields
    epmat::AbstractEPMat
    alpha::T
    beta_vec::Vector{T}
end
