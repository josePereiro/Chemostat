
function simple_toy_MetNet()
    # rxns: gt    ferm  resp  ldh   lt   biom    atpm  # mets
    S = [   1.0  -1.0   0.0   0.0   0.0   0.0    0.0;  #  G
            0.0   2.0  18.0   0.0   0.0  -55.0  -5.0;  #  E
            0.0   2.0  -1.0  -1.0   0.0   0.0    0.0;  #  P
            0.0   0.0   0.0   1.0   1.0   0.0    0.0;  #  L
        ]
    
    mets = ["G", "E", "P", "L"]
    b =    [0.0, 0.0, 0.0, 0.0] # demand
    
    metNames = ["Glucose", "Energy", "Intermediate Product" , "Lactate"];
    
    rxns = ["gt"  ,"ferm" ,"resp" , "ldh" ,  "lt" , "biom", "atpm"];
    lb =   [0.0   , 0.0   , 0.0   ,  0.0  , -100.0,   0.0,     0.5];
    ub =   [10.0 , 100.0 , 100.0 , 100.0 ,    0.0, 100.0,    100.0];
    rxnNames = ["Glucose transport", "Fermentation", "Respiration", 
        "Lactate DH", "Biomass production rate", "atp demand"];
    
    return MetNet(S, b, lb, ub, rxns, mets = mets, 
        metNames = metNames, rxnNames = rxnNames)
end

model = simple_toy_MetNet();
