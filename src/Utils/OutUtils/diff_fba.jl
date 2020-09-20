diff_va_fba_ep(bundle::ChstatBundle, ξ, β) = va_fba(bundle, ξ) - va_ep(bundle, ξ, β)
diff_va_fba_ep(bundle::ChstatBundle, ξ, β, ider) = va_fba(bundle, ξ, ider) .- va_ep(bundle, ξ, β, ider)

diff_va_fba_hr(bundle::ChstatBundle, ξ, β) = va_fba(bundle, ξ) - va_hr(bundle, ξ, β)
diff_va_fba_hr(bundle::ChstatBundle, ξ, β, ider) = va_fba(bundle, ξ, ider) .- va_hr(bundle, ξ, β, ider)
