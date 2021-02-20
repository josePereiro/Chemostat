function update_solution!(epmodel::EPModel, src_epfield::EPFields)
    dest_epfield = epmodel.epfields
    for f in fieldnames(EPFields)
        src = getfield(src_epfield, f)
        dest = getfield(dest_epfield, f)
        copyto!(dest, src)
    end
end
    
function update_solution!(epmodel::EPModel, epout::EPout)
    src_epfield = epout.sol
    T = typeof(src_epfield)
    !(T<:EPFields) && error(string("epout solution type is not valid, get '", T, "' expected ", EPFields))

    return update_solution!(epmodel, src_epfield)
end

update_solution!(epmodel_dest::EPModel, epmodel_source::EPModel) = 
    update_solution!(epmodel_dest, epmodel_source.epfields)
    