# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

## Record the state of HR algorithm
struct HRAlg{T<:Real}
    S::AbstractMatrix{T} # Stoichiometric matrix M x N sparse
    n::Integer # number of reactions
    b::AbstractVector{T} # right hand side of equation  S ν = b 
    lb::AbstractVector{T} # fluxes lower bound M elements vector
    ub::AbstractVector{T} # fluxes upper bound M elements vector
    x::AbstractVector{T} # current position
    base::AbstractMatrix{T} # nullspace(Matrix(S))
    k::Integer # number of bases
    v::AbstractVector{T} # direction holder
    
    function HRAlg(S, b, lb, ub)
        # initialization
        n = size(S, 2)
        x = warmup(S,b,lb,ub);
        base = nullspace(Matrix(S))
        v = zeros(n)
        new{eltype(S)}(S, n, b, lb, ub, x, base, size(base,2), v)
    end
end

HRAlg(model::MetNet) = HRAlg(model.S, model.b, model.lb, model.ub)