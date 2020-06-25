module LP

# TODO MathProgBase is deprecated, update the code!!!
# Prossibly make pull request for MetabolicEP too
import MathProgBase.HighLevelInterface: linprog
import Clp: ClpSolver
import ..Utils: rxnindex, FBAout, MetNet, IDER_TYPE

include("fba.jl")
include("fva.jl")


end

