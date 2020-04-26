const IDER_TYPE = Union{AbstractString, Integer}

metindex(metnet::MetNet, ider::Int) = 
    0 < ider <= length(metnet.mets) ? ider : 
        error("Index ($ider) out of range, mets [1:$(length(metnet.mets))]")
metindex(metnet::MetNet, ider::AbstractString) = 
    (indx = findfirst(isequal(ider), metnet.mets); 
    isnothing(indx) ? error("ider ($ider) not found!!!") : indx)

rxnindex(metnet::MetNet, ider::Int) = 
    0 < ider <= length(metnet.rxns) ? ider : 
        error("Index ($ider) out of range, rxns [1:$(length(metnet.rxns))]")
rxnindex(metnet::MetNet, ider::AbstractString) = 
    (indx = findfirst(isequal(ider), metnet.rxns); 
    isnothing(indx) ? error("ider ($ider) not found!!!") : indx)
