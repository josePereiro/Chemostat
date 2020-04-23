size(hrout::HRout) = (hrout.W, hrout.N)
size(hrout::HRout, i::Integer) = i == 1 ? hrout.W : i == 2 ? hrout.N : 1


