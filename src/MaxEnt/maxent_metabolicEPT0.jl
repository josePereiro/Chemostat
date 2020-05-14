# TODO add obj_ider::IDER_TYPE
function maxent_metabolicEP(model::MetNet, epout::EPout, 
    epmat::EPMatT0, obj_ider::IDER_TYPE, β; 
    maxvar=1e50,   # maximum numerical variance
    minvar=1e-50, # minimum numerical variance
    verbose = true)

M = size(epmat.G, 1)
N = length(epout.av)

G = epmat.G
Y = epmat.Y
a = epout.sol.a
b = epout.sol.b

idxy = 1:M # dependent variables
idxw = M+1:N # independent variables

ay, aw = view(a, idxy), view(a, idxw); # dep an ind prior mean (epfields)
by, bw = view(b, idxy), view(b, idxw); # dep an ind prior variance (epfields)

Σw = inv(Diagonal(1.0 ./ bw) + G' * Diagonal( 1.0 ./ by ) * G);
Σy = G*Σw*G';
vw = Σw * (aw ./ bw - G'*(ay ./ by));
vy = -G*vw .+ Y;

model_obj_idx = rxnindex(model, obj_ider)
epmat_obj_idx = findfirst(isequal(model_obj_idx), epmat.idx[1:M])

isnothing(epmat_obj_idx) && error("obj_ider '$(obj_ider)' not found between dependent fluxes. Try move it to the first index!!!")

βs_ = zeros(M);
βs_[epmat_obj_idx] = β

wy = vy + Σy*βs_
ww = pinv(-G)*(wy - Y)

# μ and σ
σy = []
μy = []

for i in 1:M # (dep (y))
    Σ_, a_, b_, v_, lb_, ub_ = Σy[i,i], ay[i], by[i], wy[i], epmat.lb[i], epmat.ub[i]
    σ = clamp(inv(1.0/Σ_ - 1.0/b_), minvar, maxvar)
    μ = Σ_ != b_ ? σ * (v_/Σ_ - a_/b_) : 0.5 * (ub_ + lb_)
    push!(σy,  σ)
    push!(μy,  μ)
end

# μ and σ
σw = []
μw = []

for i in 1:N - M # (ind (y))
    Σ_, a_, b_, v_, lb_, ub_ = Σw[i,i], aw[i], bw[i], ww[i], epmat.lb[M + i], epmat.ub[M + i]
    σ = clamp(inv(1.0/Σ_ - 1.0/b_), minvar, maxvar)
    μ = Σ_ != b_ ? σ * (v_/Σ_ - a_/b_) : 0.5 * (ub_ + lb_)
    push!(σw,  σ)
    push!(μw,  μ)
end

# epmat.idx map to the original permutation
maxflux = max(maximum(abs.(model.lb)), maximum(abs.(model.ub)))
μ = [μy; μw][epmat.idx] * maxflux;
σ = [σy; σw][epmat.idx] * maxflux^2;
    
tns = Truncated.(Normal.(μ, sqrt.(σ)), model.lb, model.ub)

return typeof(epout)(μ, σ, mean.(tns), var.(tns), epout.sol, epout.status)


end