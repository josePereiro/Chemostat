module FBA

# TODO MathProgBase is deprecated, update the code!!!
# Prossibly make pull request for MetabolicEP too
import MathProgBase.HighLevelInterface: linprog
import Clp: ClpSolver
import ..Utils: rxnindex, FBAout
import MetabolicEP: MetNet

include("fba_.jl")

end
