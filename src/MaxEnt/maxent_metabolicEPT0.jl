# TODO add obj_ider::IDER_TYPE
function maxent_metabolicEP(model::MetNet, epout::EPout, 
    epmat::EPMatT0, obj_ider::IDER_TYPE, β; 
    maxvar=1e50,   # maximum numerical variance
    minvar=1e-50, # minimum numerical variance
    verbose = true)

M = size(epmat.G, 1)
N = length(epout.av)

#scale factor
maxflux = max(maximum(abs.(model.lb)), maximum(abs.(model.ub)))

idxy = 1:M # dependent variables
idxw = M+1:N # independent variables

a = epout.sol.a .* maxflux
b = epout.sol.b .* maxflux^2
    
# TODO pull request asking for updating Σy in epmat
vw, vy, G, Y, lb, ub = epmat.vw, epmat.vy, epmat.G, epmat.Y, epmat.lb, epmat.ub

vw = vw .* maxflux
vy = vy .* maxflux
Y = Y .* maxflux
lb = lb .* maxflux
ub = ub .* maxflux

Σw = inv(Diagonal(1.0 ./ b[idxw]) + G' * Diagonal( 1.0 ./ b[idxy]) * G);
Σy = (G*Σw)*G'

fullΣ = zeros(N,N)
fullΣ[1:M, 1:M] .= Σy;
fullΣ[M + 1:N,M + 1:N] .= Σw;


obj_idx = rxnindex(model, obj_ider) # index in original model
epmat_obj_idx = findfirst(isequal(obj_idx), epmat.idx);

@show epmat_obj_idx, model.rxns[obj_idx], sum(abs.(fullΣ[:, epmat_obj_idx]))

βs_ = zeros(N);
βs_[epmat_obj_idx] = β

w = [vy; vw] + fullΣ * βs_;

# From Annas code
# Σnn = []
# wn = []
# for i in 1:N # (ind (y))
#     Σ_, a_, b_, w_, lb_, ub_ = fullΣ[i,i], a[i], b[i], w[i], lb[i], ub[i]
#     σ = clamp(inv(1.0/Σ_ - 1.0/b_), minvar, maxvar)
#     μ = Σ_ != b_ ? σ * (w_/Σ_ - a_/b_) : 0.5 * (ub_ + lb_)
#     push!(Σnn,  σ)
#     push!(wn,  μ)
# end

# From Cossios code
Σnn = b .* diag(fullΣ) ./ (b .- diag(fullΣ)) # variances of the non-truncated marginals
wn = (b .* w - diag(fullΣ) .* a) ./ (b .- diag(fullΣ)) # mean of the non-truncated marginals
    
tns = Truncated.(Normal.(wn, sqrt.(Σnn)), lb, ub)

# scale factor

# return typeof(epout)(wn[epmat.idx] .* maxflux, 
#                     Σnn[epmat.idx] .* maxflux^2, 
#                     mean.(tns)[epmat.idx] .* maxflux, 
#                     var.(tns)[epmat.idx] .* maxflux^2, 
#                     epout.sol, epout.status)
return typeof(epout)(wn[epmat.idx], 
                    Σnn[epmat.idx], 
                    mean.(tns)[epmat.idx], 
                    var.(tns)[epmat.idx], 
                    epout.sol, epout.status)


end