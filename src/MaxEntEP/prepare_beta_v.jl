# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function prepare_βv(epmat::EPMat{T}, βv::AbstractVector{T}) where T
    N = length(epmat.v)
    isempty(βv) && return spzeros(T, N)
    length(βv) != N && error("βv lenght != $N")
    return sparse(βv)
end

function prepare_βv(epmat::EPMatT0{T}, βv::AbstractVector{T}) where T
        M = length(epmat.vd) # Number of dependent variables
        N = M + length(epmat.vi)

        isempty(βv) && return spzeros(T, N - M) # same length as independent variables
        length(βv) != N && error("βv lenght != $N")
        return sparse(βv[epmat.idx[(M + 1):end]])
end

