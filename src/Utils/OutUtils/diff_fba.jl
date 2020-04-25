diff_va_fba_ep(boundle::ChstatBoundle, ξ, β) = va_fba(boundle, ξ) - va_ep(boundle, ξ, β)
diff_va_fba_ep(boundle::ChstatBoundle, ξ, β, ider) = va_fba(boundle, ξ, ider) .- va_ep(boundle, ξ, β, ider)

diff_va_fba_hr(boundle::ChstatBoundle, ξ, β) = va_fba(boundle, ξ) - va_hr(boundle, ξ, β)
diff_va_fba_hr(boundle::ChstatBoundle, ξ, β, ider) = va_fba(boundle, ξ, ider) .- va_hr(boundle, ξ, β, ider)
