module Helpers
import ..Constants: calc_fH

# Calculate pH on all scales; used in calculator.jl
"""
Calculates pH on all scales when provided with an initial pH.
"""
function calc_pH_scale(pHtot, pHfree, pHsws, pHNBS, ST, FT, temp_c, sal, Ks)
    npH = count(!isnothing, (pHtot, pHfree, pHsws, pHNBS))
    if npH != 1
        return ()
    end

    sws_to_tot_fac = (1 + ST / Ks.KS) / (1 + ST / Ks.KS + FT / Ks.KF)
    
    # Exact CO2SYS conversion factors
    offset_sws  = -log10(sws_to_tot_fac)
    offset_free = -log10(1 + ST / Ks.KS)
    # fH = a(H)_NBS / m(H)_SWS, so pHNBS = pHsws - log10(fH). pHsws is itself
    # pHtot - offset_sws, which is why the SWS-to-total term *subtracts* here. Adding it puts
    # pHNBS out by twice the SWS/total offset — 0.02 pH at S=35 — and disagrees with the
    # carbon path's `pHsws - log10(fH)`.
    offset_nbs = log10(calc_fH(temp_c, sal)) + offset_sws

    # 1. Base everything off pHtot
    local_pHtot = 0.0
    if !isnothing(pHtot)
        local_pHtot = pHtot
    elseif !isnothing(pHsws)
        local_pHtot = pHsws + offset_sws
    elseif !isnothing(pHfree)
        local_pHtot = pHfree + offset_free
    elseif !isnothing(pHNBS)
        local_pHtot = pHNBS + offset_nbs
    end

    # 2. Return all scales (reversing the offset)
    return (
        pHtot  = local_pHtot,
        pHfree = local_pHtot - offset_free,
        pHsws  = local_pHtot - offset_sws,
        pHNBS  = local_pHtot - offset_nbs
    )
end

export calc_pH_scale

end # module