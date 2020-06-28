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
    for i in non_ignored

        show_progress = verbose && (i == 1 || i % upfrec == 0 || i == n)
        show_progress && (print("fva_processing [$i / $n]        \r"); flush(stdout))

        fvalb[i], fvaub[i] = fva(S, b, fvalb, fvaub, i) .|> first

    end

    verbose && (println("fva_processing done!!!", " "^50); flush(stdout))

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