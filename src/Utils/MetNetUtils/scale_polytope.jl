function scale_polytope!(model::MetNet, factor::Real)

    for (rxni, rxn) in enumerate(model.rxns)
        lb, ub = bounds(model, rxni)
        lb -= abs(lb) * factor
        ub += abs(ub) * factor
        bounds!(model, rxni, lb, ub)
    end
    return model
end