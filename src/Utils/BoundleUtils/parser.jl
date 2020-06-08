parse_ξ(boundle, ξ::Real) = 
    ξ in boundle.ξs ? ξ : ξ isa Int ? boundle.ξs[ξ] : error("ξ: $ξ not found!!!")
parse_β(boundle, β::Real) = 
    β in boundle.βs ? β : β isa Int ? boundle.βs[β] : error("β: $β not found!!!")
