eachsample(hrout::HRout) = eachrow(hrout.hrsamples)
eachflx(hrout::HRout) = eachcol(hrout.hrsamples)
