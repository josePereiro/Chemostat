
function prepare_βv(N, βv::AbstractVector{T}) where T
    isempty(βv) && return spzeros(T, N)
    length(βv) != N && error("βv lenght != $N")
    return sparse(βv)
end
