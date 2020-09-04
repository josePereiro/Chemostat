# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function scaleepfield!(epfield, scalefact, vals...)
    @extract epfield : μ s av va 
    rmul!(μ,scalefact)
    rmul!(s,scalefact^2)
    rmul!(av,scalefact)
    rmul!(va,scalefact^2)

    for val in vals
        rmul!(val, scalefact)
    end
    
end