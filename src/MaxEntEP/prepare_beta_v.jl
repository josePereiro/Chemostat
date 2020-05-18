# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function prepare_βv(epmat::EPMat, βv::Vector)
    N = length(epmat.v)
    isempty(βv) && return zeros(N)
    length(βv) != N && error("βv lenght != $N")
    return copy(βv)
end

function prepare_βv(epmat::EPMatT0, βv::Vector)
        M = length(epmat.vy) # Number of dependent variables
        N = length(epmat.vy) + length(epmat.vw)

        isempty(βv) && return zeros(N - M) # same length as independent variables
        length(βv) != N && error("βv lenght != $N")
        return βv[epmat.idx[M + 1:end]]
end

