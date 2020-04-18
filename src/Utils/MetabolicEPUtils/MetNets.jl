MetNet(S, b, lb, ub) = MetNet(S, b, lb, ub, ["RXN $i" for i in 1:size(S,2)])
function MetNet(S, b, lb, ub, rxns)
    M,N = size(S);
    return MetNet(
        N, # number of fluxes
        M, # number of metabolites
        S, # Stoichiometric matrix M x N sparse
        b, # right hand side of equation  S ν = b 
        zeros(Float64,N), # reaction index of biomass
        lb, # fluxes lower bound N elements vector
        ub, # fluxes upper bound N elements vector
        ["NA"], # gene names
        SparseMatrixCSC{Float64,Int64}(zeros(1,1)), # rxnGeneMat
        ["NA"], # gene-reaction rule N elements vector of strings (and / or allowed)
        ["M$i" for i in 1:M], # metabolites short-name M elements 
        ["r$i" for i in 1:N], # reactions short-name N elements
        ["MET $i" for i in 1:M], # metabolites long-names M elements
        ["X" for i in 1:M], # reactions long-names N elements
        rxns, # rxn names       
        lb .!= ub, # reversibility of reactions N elements
        ["IN"]) # cellular component of fluxes N elements
end