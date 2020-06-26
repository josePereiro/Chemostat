function lineal_S(M = 3)
    N = M + 1
    # S
    S = zeros(M,N)
    for i in 1:M
        S[i,i] = 1.0
        S[i,i+1] = -1.0
    end
    return SparseMatrixCSC{Float64,Int64}(S)
end

function simple_lineal_MetNet(M)
    S = lineal_S(M)
    M,N = size(S)
    return MetNet(S, zeros(Float64, M), zeros(Float64,N), ones(Float64, N))
end

