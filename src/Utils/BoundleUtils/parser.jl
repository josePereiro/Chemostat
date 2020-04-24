parse_ξ(chdata, ξ::Real) = ξ in chdata.ξs ? ξ : ξ isa Int ? chdata.ξs[ξ] : error("ξ: $ξ not found!!!")
parse_β(chdata, β::Real) = β in chdata.βs ? β : β isa Int ? chdata.βs[β] : error("β: $β not found!!!")
