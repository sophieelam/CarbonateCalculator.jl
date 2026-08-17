module Helpers
using Printf
using ForwardDiff, LinearAlgebra

"""
    calc_fH(temp_c, sal)

Activity coefficient ratio linking the NBS and seawater pH scales, `a(H)_NBS / m(H)_SWS`.

Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80), via CO2SYS.

Temperature is in Celsius, as everywhere else in the package; the polynomial itself is
fitted in Kelvin and the conversion happens here, so there is one place to get it wrong
rather than one per call site. Expect ~0.71 for seawater.

This is the only definition — `Constants` reaches it through `using ..Helpers`.
"""
function calc_fH(temp_c, sal)
    temp_k = temp_c + 273.15
    a, b, c, d = (1.2948, -2.036e-3, 4.607e-4, -1.475e-6)
    return a + b * temp_k + (c + d * temp_k) * sal^2
end


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
    # fH = a(H)_NBS / m(H)_SWS, so pHNBS = pHsws - log10(fH); since pHsws is itself
    # pHtot - offset_sws, the SWS-to-total term subtracts here. It used to be added, which
    # put pHNBS out by twice the SWS/total offset (0.02 pH at S=35) and disagreed with the
    # carbon path's `final_pHsws - log10(fH_val)`.
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


"""
Function to format outputs of boron/boron isotope calculator results
"""
function print_boron_results(ps)
    println("=========================================")
    println("        BORON SYSTEM TEST RESULTS        ")
    println("=========================================")
    println("Inputs:")
    println("  pH (Total)   : ", round(ps.pHtot, digits=3))
    println("  Total Boron  : ", round(ps.BT, digits=2))
    println("  δ11B Total   : ", round(ps.δBT, digits=2))
    println("-----------------------------------------")
    println("Calculated Concentrations:")
    println("  Boric Acid (BOH₃) : ", round(ps.BOH₃, digits=2), " µmol/kg")
    println("  Borate (BOH₄)     : ", round(ps.BOH₄, digits=2), " µmol/kg")
    println("-----------------------------------------")
    println("Calculated Isotopes:")
    println("  δ11B of BOH₃ : ", round(ps.δBOH₃, digits=2), " ‰")
    println("  δ11B of BOH₄ : ", round(ps.δBOH₄, digits=2), " ‰")
    println("=========================================")
end


"""
Function to format outputs of whole system calculator results
"""
function print_system_results(results::NamedTuple)
    println("\n==============================================")
    println("          WHOLE SYSTEM TEST RESULTS           ")
    println("==============================================")
    println(rpad("Parameter", 22), " | ", "Value")
    println("-----------------------|----------------------")
    
    for (k, v) in pairs(results)
        # We skip printing massive nested dictionaries/objects if any snuck through
        if v isa AbstractFloat
            # Floats get rounded to 4 decimal places for clean reading
            @printf("%-22s | %10.4f\n", string(k), v)
        elseif !isnothing(v) && !(v isa Dict) && !(v isa NamedTuple)
            @printf("%-22s | %10s\n", string(k), string(v))
        elseif isnothing(v)
            @printf("%-22s | %10s\n", string(k), "missing")
        end
    end
    println("==============================================\n")
end

export calc_fH, calc_pH_scale

end # module