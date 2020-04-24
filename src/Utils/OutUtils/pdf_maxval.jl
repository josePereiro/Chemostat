function pdf_maxval(model::MetNet, out::Union{FBAout, EPout}, ider)
    tN = marginal(model, out, ider)
    pdf_maxval_ = max(pdf(tN, mean(tN)), pdf(tN, lb(model, ider)), pdf(tN, ub(model, ider)))
    return -Inf < pdf_maxval_ < Inf ? pdf_maxval_ : -1
end

function pdf_maxval(model::MetNet, out::HRout, ider)
    hist = hists(model, out, ider)
    hist = normalize(hist, mode = :pdf)
    return maximum(hist.weights)
end

function pdf_maxval(model::MetNet, outs, ider)
    pdf_maxval_ = maximum([pdf_maxval(model, out, ider) for out in outs])
    return -Inf < pdf_maxval_ < Inf ? pdf_maxval_ : -1
end