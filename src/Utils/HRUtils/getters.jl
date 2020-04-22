flx_samples(model::MetNet, hrsamples::Matrix, ider) = hrsamples[:, rxnindex(model, ider)]

va(model::MetNet, hrsamples::Matrix, ider) = var(flx_samples(model, hrsamples, ider))
av(model::MetNet, hrsamples::Matrix, ider) = mean(flx_samples(model, hrsamples, ider))