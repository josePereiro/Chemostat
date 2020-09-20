parse_ξ(bundle, ξ::Real) = 
    ξ in bundle.ξs ? ξ : ξ isa Int ? bundle.ξs[ξ] : error("ξ: $ξ not found!!!")
parse_β(bundle, β::Real) = 
    β in bundle.βs ? β : β isa Int ? bundle.βs[β] : error("β: $β not found!!!")
