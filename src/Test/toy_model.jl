function simple_toy_MetNet()
    net = Dict()
    net[:S] = 
    # rxns: gt    ferm  resp  ldh   lt   biom    atpm  # mets
        [   1.0  -1.0   0.0   0.0   0.0   0.0    0.0;  #  G
            0.0   2.0  18.0   0.0   0.0  -55.0  -5.0;  #  E
            0.0   2.0  -1.0  -1.0   0.0   0.0    0.0;  #  P
            0.0   0.0   0.0   1.0   1.0   0.0    0.0;  #  L
        ]
    
    net[:mets] = ["G", "E", "P", "L"]
    net[:b] =    [0.0, 0.0, 0.0, 0.0] # demand
    
    net[:metNames] = ["Glucose", "Energy", "Intermediate Product" , "Lactate"];
    
    net[:rxns] = ["gt"  ,"ferm" ,"resp" , "ldh" ,  "lt" , "biom", "atpm"];
    net[:lb] =   [0.0   , 0.0   , 0.0   ,  0.0  , -100.0,   0.0,     0.5];
    net[:ub] =   [10.0 , 100.0 , 100.0 , 100.0 ,    0.0, 100.0,    100.0];
    net[:c] =    [ 0.0 ,   0.0 ,   0.0 ,   0.0 ,    0.0,   0.0,      0.0];
    net[:rxnNames] = ["Glucose transport", "Fermentation", "Respiration", 
        "Lactate DH", "Lactate transport", "Biomass production rate", "atp demand"];
    
    return MetNet(;net...)
end

function  toy_model(;resp_cost::Real = -0.1, E_demand::Real = 5.0)
    @assert resp_cost <= 0.0

    model = simple_toy_MetNet()

    # Atpm
    lb!(model, "atpm", E_demand)

    # cost
    # Add a new metabolite simulating the cost penalazing 
    # reaction fluxes. 
    # A new balance equations is then added:
    #     Σ(rᵢ*costᵢ) + tot_cost = 0
    # Because the cost coefficients (costᵢ) < 0, the system must allocate 
    # the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
    # are usually bounded [0.0, 1.0]

    # We will add 1 met and 1 rxn
    M, N = size(model)
    model = expanded_model(model, M + 1, N + 1)
    # Add cost met
    set_met!(model, M + 1, Met("cost"; rxns = ["resp"], S = [resp_cost]))
    set_rxn!(model, N + 1, Rxn("tot_cost"; rxns = ["cost"], S = [1.0], lb = 0.0, ub = 1.0))
    return fva_preprocess(model)
end