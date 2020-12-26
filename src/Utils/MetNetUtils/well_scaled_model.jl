"""
    returns the min and max absolute non zero value
"""
function nzabs_range(col::AbstractArray{<:Number})
    nzabs = filter(!iszero, col) .|> abs
    isempty(nzabs) && return (zero(eltype(nzabs)), zero(eltype(nzabs)))
    return (minimum(nzabs), maximum(nzabs))
end


"""
    returns (a, n) so that s ≈ a^n. a ∈ (0, b]
"""
function factorize(s::Real, b::Float64 = 100.0)::Tuple{Float64,Int64}
    iszero(s) && return (0.0, 1)
    n = ceil(Int, abs(log(b, s)))
    a = s^(1/n)
    return (a, n)
end
factorize(s::Real, b::Real) = factorize(s, float(b))

"""
    Compute the future size of the scaled model
"""
compute_well_scaled_size(model::MetNet, b = 100.0) = 
    size(model) .+ sum(max.(last.(factorize.(abs.(model.S[model.S .!= 0]), b)) .- 1, 0))

const LIFT_TAG = "LIFT"
function well_scaled_model(model::MetNet, b::Real = 100.0; 
        lift_bound::Real = 1e50, verbose = true)

    ## Inferring final size
    eM, eN = compute_well_scaled_size(model, b)
    exp_model = expanded_model(model, eM, eN)

    verbose && (prog = Progress(size(model, 1); desc = "Scaling   ", dt = 0.5))
    for (rxni, rxn) in model.rxns |> enumerate

        metis = rxn_mets(model, rxni)
        mets = model.mets[metis]
        for (meti, met) in zip(metis, mets)

            s = model.S[meti, rxni]
            sign_s = sign(s)
            news, n = factorize(abs(s), b)

            if n > 1 
                # We need to lift
                # see (https://doi.org/10.1186/1471-2105-14-240)

                # orig rxn
                exp_model.S[meti, rxni] = 0.0

                # lift met (n - 1)
                lift_met_num = n - 1
                lift_met = string(LIFT_TAG, lift_met_num, "_", met)
                lift_meti = findempty(exp_model, :mets)
                exp_model.mets[lift_meti] = lift_met
                exp_model.S[lift_meti, rxni] = sign_s * news

                # lift rxn chain 
                lift_rxn_num = 1
                for i in 1:(n - 2)
                    # new lift rxn
                    lift_rxn = string(rxn, "_", LIFT_TAG, lift_rxn_num)
                    lift_rxni = findempty(exp_model, :rxns)
                    exp_model.rxns[lift_rxni] = lift_rxn
                    exp_model.lb[lift_rxni], exp_model.ub[lift_rxni] = -lift_bound, lift_bound
                    
                    # current lift met
                    exp_model.S[lift_meti, lift_rxni] = -sign_s * 1.0
                    lift_met_num -= 1

                    # new lift met
                    lift_met = string(LIFT_TAG, lift_met_num, "_", met)
                    lift_meti = findempty(exp_model, :mets)
                    exp_model.mets[lift_meti] = lift_met
                    exp_model.S[lift_meti, lift_rxni] = sign_s * news

                    lift_rxn_num += 1
                end

                # last lift rxn (n - 1)
                lift_rxn = string(rxn, "_", LIFT_TAG, lift_rxn_num)
                lift_rxni = findempty(exp_model, :rxns)
                exp_model.rxns[lift_rxni] = lift_rxn
                
                # s * orig_met --> lift_met1
                exp_model.S[meti, lift_rxni] = sign_s * news # same sign
                exp_model.S[lift_meti, lift_rxni] = -sign_s * 1.0
                exp_model.lb[lift_rxni], exp_model.ub[lift_rxni] = -lift_bound, lift_bound

            end # if n > 1 

        end # zip(metis, mets)

        verbose && update!(prog, rxni)
    end
    verbose && finish!(prog)

    return exp_model
end
