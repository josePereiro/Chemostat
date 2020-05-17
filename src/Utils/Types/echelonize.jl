# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function idxlicols(X; tol::Float64=1e-10)

    sum(abs2,X) == 0 && (return(Array{Int,1}(), Array{Int,2}() ))
    Q,R,E = qr(X, Val(true))
    diagr = abs.(diag(R))
    r = findall(diagr .>= tol*diagr[1])[end]
    idx = sort(E[1:r])
    return idx
end
    
function echelonize(X::T,v; eps::Real=1e-10) where T <:DenseArray
    M,N = size(X)

    idxrow = idxlicols(Matrix(X'))
    Mred = length(idxrow)

    idxind = idxlicols(X)
    idxdep = setdiff(1:N,idxind)
    newidx = vcat(idxind,idxdep)
    Tv = @view X[idxrow,newidx]
    iTv = inv(Tv[1:Mred, 1:Mred])
    res = iTv * Tv
    # trimming zeros
    for i in eachindex(res)
        abs(res[i]) < eps && (res[i] = zero(res[i]))
    end
    bnew = iTv * v[idxrow]
    # trimming ones
    for i in 1:Mred
        abs(1.0 - res[i,i]) < eps && (res[i,i] = one(res[i,i]))
    end

    return idxdep, newidx, res, bnew

end