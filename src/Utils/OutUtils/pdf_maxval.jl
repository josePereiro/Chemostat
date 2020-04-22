function pdf_maxval(model::MetNet, out::Union{FBAout, EPout}, ider)
    ider = rxnindex(model, ider)
    D = marginal(model, out, ider)
    return pdf(D, mean(D))
end

function pdf_maxval(model::MetNet, out::HRout, ider)
    ider = rxnindex(model, ider)
    hist = hists(model, out, ider)
    hist = normalize(hist, mode = :pdf)
    return maximum(hist.weights)
end

pdf_maxval(model::MetNet, outs, ider) = 
    maximum([pdf_maxval(model, out, ider) for out in outs])