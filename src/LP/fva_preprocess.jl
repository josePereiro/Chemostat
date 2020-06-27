# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
function fva_preprocess(S,b,lb,ub,rxns; 
        verbose = false, eps = 0.0, 
        ignored = [], # rxns skip preprocess
        protected = [], # rxns skip blocking
        upfrec = 10)
        
    b, lb, ub = copy(b), copy(lb), copy(ub)
    n = length(lb)
    ei = zeros(n)
    blocked = falses(n)
    upfrec = floor(Int, n/upfrec)
    for i=1:n

        i in ignored && continue
        
        show_progress = verbose && (i == 1 || i % upfrec == 0 || i == n)
        show_progress && (print("fva_processing [$i / $n]        \r"); flush(stdout))

        lb_, ub_ = fva(S, b, lb, ub, i) .|> first

        if lb_ == ub_ # blocking
            i in protected && continue
            blocked[i] = true
            lb[i], ub[i] = lb_ - eps, ub_ + eps
        else
            lb[i] = lb_
            ub[i] = ub_
        end

    end
    verbose && (println("fva_processing done!!!                                           "); flush(stdout))
    idxb = findall(blocked)

    verbose && println("$(sum(blocked)) blocked fluxes")

    unblocked = (lb .< ub)

    return S[:,unblocked], b - S*(blocked .* lb), lb[unblocked], ub[unblocked], rxns[unblocked], idxb
end


function fva_preprocess(metnet::MetNet; 
            verbose = false, eps = 0.0, return_blocked = false,
            ignored = [], # rxns skip preprocess
            protected = [] # rxns skip bloking

        )
    ignored = map((r) -> rxnindex(metnet, r), ignored)
    protected = map((r) -> rxnindex(metnet, r), protected)

    S_, b_, lb_, ub_, rxns_, idxb_ = fva_preprocess(metnet.S, metnet.b, metnet.lb, 
        metnet.ub, metnet.rxns; verbose = verbose, eps = eps, ignored = ignored);
    
    metnet = MetNet(metnet; S = S_, b = b_, lb = lb_, ub = ub_, rxns = rxns_)
    return return_blocked ? (metnet, idxb_) : metnet
end