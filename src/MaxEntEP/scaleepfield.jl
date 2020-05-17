# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function scaleepfield!(epfield,lb, ub,Y,scalefact)
    @extract epfield : μ s av va
    rmul!(μ,scalefact)
    rmul!(s,scalefact^2)
    rmul!(av,scalefact)
    rmul!(va,scalefact^2)
    rmul!(ub,scalefact)
    rmul!(lb,scalefact)
    rmul!(Y,scalefact)
end