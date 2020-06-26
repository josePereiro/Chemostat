function make_intake_info_compatible!(model::MetNet)
    for (intake, info) in model.intake_info
        if !(intake in model.rxns)
            delete!(model.intake_info, intake)
        end
    end
    return model
end