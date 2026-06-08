module Helpers
using Printf
using ForwardDiff, LinearAlgebra

# Calculates fH; used in carbon.jl.
"""
Taken from CO2SYS
Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80)
Note: Temperature MUST be in Celsius!
"""
function calc_fH(T_K, Salinity)
    a, b, c, d = (1.2948, -2.036e-3, 4.607e-4, -1.475e-6)
    return a + b * T_K + (c + d * T_K) * Salinity^2
end


# Calculate pH on all scales; used in calculator.jl
"""
Calculates pH on all scales when provided with an initial pH.
"""
function calc_pH_scale(pHtot, pHfree, pHsws, pHNBS, ST, FT, TempC, S, Ks)
    npH = count(!isnothing, (pHtot, pHfree, pHsws, pHNBS))
    if npH != 1
        return ()
    end

    sws_to_tot_fac = (1 + ST / Ks.KS) / (1 + ST / Ks.KS + FT / Ks.KF)
    
    # Exact CO2SYS conversion factors
    offset_sws  = -log10(sws_to_tot_fac)
    offset_free = -log10(1 + ST / Ks.KS)
    offset_nbs = log10(calc_fH(TempC, S)) + log10(sws_to_tot_fac)

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
Calculates total sulfur [mol/kg-sw]
Taken from CO2SYS (From Dickson et al., 2007, Table 2)
Note: Salinity / 1.80655 = Chlorinity
"""
function calc_ST(S)
    return (0.14 * S / 1.80655 / 90.062)
end 


"""
Calculates total fluorine [mol/kg-sw]
Taken from CO2SYS (From Dickson et al., 2007, Table 2)
Note: Salinity / 1.80655 = Chlorinity
"""
function calc_FT(S)
    return (6.7e-5 * S / 1.80655 / 18.9984)
end


"""
Calculates total boron [mol/kg-sw]
Taken from CO2SYS (Uppstrom, L., Deep-Sea Research 21:161-162, 1974)
this is 0.0004157 * Sal/35. = 0.0000119 * Sal
TB(FF) = (0.000232 / 10.811) * (Sal / 1.80655) in mol/kg-SW
"""
function calc_BT(S)
    return (0.000416 / 35.0) * S
end


"""
Function to format outputs of calculator results
"""
function print_carbon_results(ps)
    println("=========================================================")
    println("                CARBON SYSTEM CALCULATOR                 ")
    println("=========================================================")
    @printf("%-20s | %-15s | %-15s\n", "Parameter", "Input Cond.", "Output Cond.")
    println("---------------------------------------------------------")
    
    # Core Parameters
    @printf("%-20s | %-15.2f | %-15.2f\n", "Temperature (°C)", ps.T_in, ps.T_out)
    @printf("%-20s | %-15.2f | %-15.2f\n", "Salinity", ps.S_in, ps.S_out)
    @printf("%-20s | %-15.2f | %-15.2f\n", "Pressure", ps.P_in, ps.P_out)
    println("---------------------------------------------------------")
    
    # Carbon Species (using the unit from the bundle)
    unit = ps.unit
    @printf("%-20s | %-15.2f | %-15.2f\n", "DIC ($unit)", ps.DIC, ps.DIC_in)
    @printf("%-20s | %-15.2f | %-15.2f\n", "TA ($unit)", ps.TA, ps.TA_in)
    @printf("%-20s | %-15.2f | %-15.2f\n", "pCO₂ (uatm)", ps.pCO₂, ps.pCO₂_in)
    @printf("%-20s | %-15.2f | %-15.2f\n", "HCO₃ ($unit)", ps.HCO₃, ps.HCO₃_in)
    @printf("%-20s | %-15.2f | %-15.2f\n", "CO₃ ($unit)", ps.CO₃, ps.CO₃_in)
    println("---------------------------------------------------------")
    
    # pH and Saturation
    @printf("%-20s | %-15.4f | %-15.4f\n", "pH (Total Scale)", ps.pHtot, ps.pHtot_in)
    @printf("%-20s | %-15.4f | %-15.4f\n", "Ω Aragonite", ps.ΩA, ps.ΩA_in)
    @printf("%-20s | %-15.4f | %-15.4f\n", "Ω Calcite", ps.ΩC, ps.ΩC_in)
    @printf("%-20s | %-15.4f | %-15.4f\n", "Revelle Factor", ps.revelle_factor, ps.revelle_factor_in)
    println("=========================================================")
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
export calc_fH, calc_ST, calc_FT, calc_pH_scale, calc_BT

end # module