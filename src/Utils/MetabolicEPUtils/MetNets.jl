fake_metsid(M) = ["M$i" for i in 1:M]
fake_metNames(M) = ["MET $i" for i in 1:M]
fake_metFormulas(M) = ["X" for i in 1:M]


fake_rxnsid(N) = ["r$i" for i in 1:N]
fake_rxnNames(N) = ["RXN $i" for i in 1:N]

        
function MetNet(S, b, lb, ub, rxns = fake_rxnsid(size(S,2));
    mets = fake_metsid(size(S,1)),  
    metNames =  fake_metNames(size(S,1)), 
    metFormulas = fake_metFormulas(size(S,1)), 
    rxnNames = fake_rxnNames(size(S,2)))
    
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
        mets, # metabolites short-name M elements 
        rxns, # reactions short-name N elements
        metNames, # metabolites long-names M elements
        metFormulas, # metabolites formula M elements
        rxnNames, # reactions long-name N elements
        (lb .!= 0.0) .& (ub .!= 0.0), # reversibility of reactions N elements
        ["IN"]) # cellular component of fluxes N elements
end