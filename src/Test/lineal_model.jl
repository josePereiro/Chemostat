function lineal_S(M = 3)
    N = M + 1
    # S
    S = zeros(M,N)
    for i in 1:M
        S[i,i] = 1.0
        S[i,i+1] = -1.0
    end
    return S
end

function simple_lineal_MetNet(M)
    net = Dict()
    net[:S] = lineal_S(M)
    M,N = size(net[:S])
    net[:c] = zeros(N)
    net[:b] = zeros(M)
    net[:lb] = zeros(N)
    net[:ub] = ones(N)

    return MetNet(net; reshape = false)
end

