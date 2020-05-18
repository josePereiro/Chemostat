struct FBAout{T<:Real} <: AbstractOut
    v::Vector{T} # Flux vector solution
    obj_val::T # The value of the objective
    obj_ider # The used objective identifier
    sol # The LP solution object 
end