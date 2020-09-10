# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
function fva_preprocess(S,b,lb,ub,rxns; 
        verbose = false, eps = 0.0, 
        ignored = [], # rxns skip preprocess
        protected = [], # rxns skip blocking
        upfrec = 10)
        
    fvalb, fvaub = (lb, ub) .|> copy

    m, n = size(S)
    ei = zeros(n)
    blocked = falses(n)
    upfrec = floor(Int, n/upfrec)

    
    # FVA
    _bidx = trues(n)
    _bidx[ignored] .= false
    non_ignored = findall(_bidx)
    nn = length(non_ignored) 

    prog = Progress(nn; desc = "Doing FVA  ")
    for i in non_ignored
        fvalb[i], fvaub[i] = fva(S, b, fvalb, fvaub, i; verbose = false) .|> first
        verbose && update!(prog, i)
    end

    verbose && finish!(prog)

    return del_blocked(S, b, fvalb, fvaub, rxns; 
                    eps = eps, protected = protected)
end


function fva_preprocess(metnet::MetNet; 
            verbose = false, eps = 0.0, return_blocked = false,
            ignored = [], # rxns skip preprocess
            protected = [] # rxns skip blocking
        )
    ignored = map((r) -> rxnindex(metnet, r), ignored)
    protected = map((r) -> rxnindex(metnet, r), protected)

    S_, b_, lb_, ub_, rxns_, blocked = fva_preprocess(metnet.S, metnet.b, metnet.lb, 
        metnet.ub, metnet.rxns; verbose = verbose, eps = eps, ignored = ignored);
    
    metnet = MetNet(metnet; S = S_, b = b_, lb = lb_, ub = ub_, rxns = rxns_)
    return return_blocked ? (metnet, blocked) : metnet
end