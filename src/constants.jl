"""
    Constants

Equilibrium constants, their temperature/salinity/pressure dependence, and the seawater
composition derived from salinity.

Every constant is available under several published parameterisations, registered in
[`COMPOSITION`](@ref), [`SPECIATION`](@ref) and [`CONVERSION`](@ref) and chosen with a
[`KSettings`](@ref). [`calculate_constants`](@ref) assembles a complete set for one body of
water.

Each choice is independent: a K1/K2 parameterisation sets K1 and K2 and nothing else. Where a
parameterisation was published or used alongside particular companions (GEOSECS, Peng et al.,
the freshwater constants), the registry records them and `KSettings` warns when they are not
selected, rather than switching them silently.

Constants are returned on the total pH scale, pressure-corrected after Millero (1995).
"""
module Constants
using Kgen

# Global Constants

"The gas constant in ml bar⁻¹ K⁻¹ mol⁻¹. DOEv2 gives 83.14462618; this matches CO2SYS."
const GAS_CONSTANT = 83.14472

# Modern seawater composition, taken from Kgen so there is a single source of truth for
# what "unspecified Ca/Mg" means. Ca and Mg default to `nothing` throughout the public
# API: nothing leaves Ca_method in charge of the calcium used for saturation state, and
# means "modern" for the KGen correction.
const MODERN_CALCIUM = Kgen.MODERN_CALCIUM
const MODERN_MAGNESIUM = Kgen.MODERN_MAGNESIUM

"Return the factor converting a seawater-scale constant to the total scale."
function SWStoTOT(; ST, FT, KS, KF, kwargs...)
    return (1 + ST / KS) / (1 + ST / KS + FT / KF)
end

"Return the factor converting a free-scale constant to the total scale."
function FREEtoTOT(; ST, KS, kwargs...)
    return (1 + ST/KS)
end

# 1. Roy et al., 1993: https://doi.org/10.1016/0304-4203(93)90207-5
# Salinity between 5-45 PSU, temperature between 0-45 C
# Total pH scale, artificial seawater
# NOTES FROM CO2SYS:
# Typo: in the abstract on p. 249: in the eq. for lnK1* the
# last term should have S raised to the power 1.5.
# They claim standard deviations (p. 254) of the fits as
# .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
# They also claim (p. 258) 2s precisions of .004 in pK1 and
# .006 in pK2. These are consistent, but Andrew Dickson
# (personal communication) obtained an rms deviation of about
# .004 in pK1 and .003 in pK2. This would be a 2s precision
# of about 2% in K1 and 1.5% in K2.
# T:  0-45  S:  5-45. Total Scale. Artificial sewater.
# This is eq. 29 on p. 254 and what they use in their abstract:
function Roy1993(; temp_c, sal, ST, FT, KS, KF, kwargs...)
    TK = temp_c + 273.15
    logTK = log(TK)
    sqrS = sqrt(sal)
    S15 = sal * sqrS

    # K1 Calculation
    lnK1 = (2.83655 - 2307.1266 / TK - 1.5529413 * logTK +
           (-0.20760841 - 4.0484 / TK) * sqrS + 0.08468345 * sal -
           0.00654208 * S15)
    
    # Originally on total scale! "* (1 - 0.001005109 * sal)" converts to mol/kg-sw
    # then, if I were to "/ SWStoTOT" I would get to SWS scale
    K1 = exp(lnK1) * (1 - 0.001005109 * sal)

    # K2 Calculation 
    lnK2 = (-9.226508 - 3351.6106 / TK - 0.2005743 * logTK +
           (-0.106901773 - 23.9722 / TK) * sqrS + 0.1130822 * sal -
           0.00846934 * S15)
    
    # Originally on total scale! "* (1 - 0.001005109 * sal)" converts to mol/kg-sw
    # then, if I were to "/ SWStoTOT" I would get to SWS scale
    K2 = exp(lnK2) * (1 - 0.001005109 * sal)

    return (; K1=K1, K2=K2)
end


# 2. Goyet & Poisson, 1989: https://doi.org/10.1016/0198-0149(89)90064-2
# Salinity between 10-50 PSU, temperature between -1-40 C
# Seawater pH scale, artificial seawater
# NOTES FROM CO2SYS:
# The 2s precision in pK1 is .011, or 2.5% in K1.
# The 2s precision in pK2 is .02, or 4.5% in K2.
# This is in Table 5 on p. 1652 and what they use in the abstract:
function GP1989(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = (812.27 / TK + 3.356 - 0.00171 * sal * log(TK) +
    0.000091 * sal^2)

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = (1450.87 /TK + 4.604 - 0.00385 * sal * log(TK) +
    0.000182 * sal^2)

    K2 = 10.0^(-pK2)

    return (; K1=K1, K2=K2)
end


# 3. Hansson, 1973: https://doi.org/10.1016/0011-7471(73)90100-9
# Salinity between 20-40 PSU, temperature between 2-35 C
# Seawater pH scale, artificial seawater
# NOTES FROM CO2SYS:
# HANSSON refit BY DICKSON AND MILLERO
# Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
# (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
# refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
# and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
# on the SWS pH scale in mol/kg-SW.
# Hansson gave his results on the Total scale (he called it
# the seawater scale) and in mol/kg-SW.
# Typo in DM on p. 1739 in Table 4: the equation for pK2*
# for Hansson should have a .000132 *S^2
# instead of a .000116 *S^2.
# The 2s precision in pK1 is .013, or 3% in K1.
# The 2s precision in pK2 is .017, or 4.1% in K2.
# This is from Table 4 on p. 1739.
function Hansson1973(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = 851.4 /TK + 3.237 - 0.0106 * sal + 0.000105 * sal^2

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = (-3885.4 / TK + 125.844 - 18.141 * log(TK) - 0.0192 * sal +
    0.000132 * sal^2)

    K2 = 10.0^(-pK2)

    return (; K1=K1, K2=K2)
end 


# 4. Dickson & Millero, 1987: https://doi.org/10.1016/0198-0149(87)90021-5
# Salinity between 20-40 PSU, temperature between 2-35 C
# Seawater pH scale, artificial seawater
# NOTES FROM CO2SYS
# (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
# refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
# on the SWS pH scale in mol/kg-SW.
# Mehrbach et al gave results on the NBS scale.
# The 2s precision in pK1 is .011, or 2.6% in K1.
# The 2s precision in pK2 is .020, or 4.6% in K2.
# Valid for salinity 20-40.
# This is in Table 4 on p. 1739.
function DM1987(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = (3670.7 / TK - 62.008 + 9.7944 * log(TK) - 0.0118 * sal +
    0.000116 * sal^2)

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = 1394.7 / TK + 4.777 - 0.0184 * sal + 0.000118 * sal^2

    K2 = 10.0^(-pK2)

    return (; K1=K1, K2=K2)
end 


# 5. Hansson & Mehrbach 1973: https://doi.org/10.3891/ACTA.CHEM.SCAND.27-0931
# Salinity between 20-40 PSU, temperature between 2-35 C
# Seawater pH scale, artificial seawater
# NOTES FROM CO2SYS:
# HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
# Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
# (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
# refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
# Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
# and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
# on the SWS pH scale in mol/kg-SW.
# Typo in DM on p. 1740 in Table 5: the second equation
# should be pK2* =, not pK1* =.
# The 2s precision in pK1 is .017, or 4% in K1.
# The 2s precision in pK2 is .026, or 6% in K2.
# Valid for salinity 20-40.
# This is in Table 5 on p. 1740.
function HM1973(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = 845 / TK + 3.248 - 0.0098 * sal + 0.000087 * sal^2

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = 1377.3 / TK + 4.824 - 0.0185 * sal + 0.000122 * sal^2

    K2 = 10.0^(-pK2)

    return (; K1=K1, K2=K2)
end


# 6/7. Mehrbach 1973: https://doi.org/10.4319/lo.1973.18.6.0897
# Salinity between 19-43 PSU, temperature between 2-35 C
# NBS pH scale converted to SWS, real seawater
# NOTES FROM CO2SYS:
# GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
# Limnology and Oceanography, 18(6):897-907, 1973.
# I.e., these are the original Mehrbach dissociation constants.
# The 2s precision in pK1 is .005, or 1.2% in K1.
# The 2s precision in pK2 is .008, or 2% in K2.
function Mehrbach1973(; temp_c, sal, fH, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = (-13.7201 + 0.031334 * TK + 3235.76 / TK + 1.3e-5 * sal *
    TK - 0.1032 * sqrt(sal))

    # On the SWS pH scale in mol/kg-sw (converted from NBS scale)
    K1 = 10.0^(-pK1) / fH

    pK2 = (5371.9645 + 1.671221 * TK + 0.22913 * sal + 18.3802 *
    log10(sal)) - 128375.28 / TK - 2194.3055 * log10(TK) -
    8.0944e-4 * sal * TK - 5617.11 * log10(sal) / TK + 2.136 *
    sal / TK

    K2 = 10.0^(-pK2) / fH

    return (; K1=K1, K2=K2)
end


# 8. Millero, 1979: https://doi.org/10.1016/0016-7037(79)90184-4
# Salinity of 0 PSU (freshwater), temperature between 0-50 C
# NOTES FROOM CO2SYS:
# K1 from refit data from Harned and Davis,
# J American Chemical Society, 65:2030-2037, 1943.
# K2 from refit data from Harned and Scholes,
# J American Chemical Society, 43:1706-1709, 1941.
# This is only to be used for Sal=0 water (note the absence of S in the
# below formulations).
# These are the thermodynamic Constants:
function Millero1979(; temp_c, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    lnK1 = 290.9097 - 14554.21 / TK - 45.0575 * log(TK)

    K1 = exp(lnK1)

    lnK2 = 207.6548 - 11843.79 / TK - 33.6485 * log(TK)

    K2 = exp(lnK2)

    return (; K1=K1, K2=K2)
end 


# 9. Cai & Wang 2003: https://doi.org/10.4319/lo.1998.43.4.0657
# Salinity between 0-49 PSU, temperature between 2-35 C
# NBS pH scale converted to SWS, artificial and real seawater
# For estuarine use
# NOTES FROM CO2SYS:
# Data used in this work is from:
# K1: Merhback (1973) for S>15, for S<15: Mook and Keone (1975)
# K2: Merhback (1973) for S>20, for S<20: Edmond and Gieskes (1970)
# Sigma of residuals between fits and above data: Â±0.015, +0.040 for K1
# and K2, respectively.
# Sal 0-40, Temp 0.2-30
# Limnol. Oceanogr. 43(4) (1998) 657-668
# On the NBS scale
# Their check values for F1 don't work out, not sure if this was correctly
# published...
# Conversion to SWS scale by division by fH is uncertain at low Sal due to
# junction potential.
function CW2003(; temp_c, sal, fH, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin
    F1 = 200.1 / TK + 0.3220
    F2 = -129.24 / TK + 1.4381

    pK1 = 3404.71 / TK + 0.032786 * TK - 14.8435 - 0.071692 * F1 * 
    sqrt(sal) + 0.0021487 * sal

    K1 = 10.0^(-pK1) / fH

    pK2 =  2902.39 / TK + 0.02379 * TK - 6.4980 - 0.3191 * F2 * 
    sqrt(sal) + 0.0198 * sal

    K2 = 10.0^(-pK2) / fH

    return (; K1=K1, K2=K2)
end 


# 10. Lueker, et al., 2000: https://doi.org/10.1016/S0304-4203(00)00022-0
# Salinity between 19-43 PSU, temperature between 2-35 C
# Total pH scale, real seawater (converted to SWS)
# NOTES FROM CO2SYS:
# This is Mehrbach's data refit after conversion to the Total scale, for
# comparison with their equilibrator work.
# Mar. Chem. 70 (2000) 105-119
# Total scale and kg-sw
function Lueker2000(; temp_c, sal, ST, FT, KS, KF, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = 3633.86 / TK - 61.2172 + 9.6777 * log(TK) - 0.011555 *
    sal + 0.0001152 * sal^2

    K1 = 10.0^(-pK1) # / SWStoTOT(ST=ST, FT=FT, KS=KS, KF=KF) # converted to SWS

    pK2 =  471.78 / TK + 25.929 - 3.16967 * log(TK) - 0.01781 *
    sal + 0.0001122 * sal^2

    K2 = 10.0^(-pK2) # / SWStoTOT(ST=ST, FT=FT, KS=KS, KF=KF) # converted to SWS

    return (; K1=K1, K2=K2) # By commenting /SWStoTOT out, on total
end 


# 11. Mojica Prieto & Millero, 2002: https://doi.org/10.1016/S0016-7037(02)00855-4
# Salinity between 5-42 PSU, temperature between 0-45 C
# pH on seawater scale, real seawater
# NOTES FROM CO2SYS:
# sigma for pK1 is reported to be 0.0056
# sigma for pK2 is reported to be 0.010
# This is from the abstract and pages 2536-2537
function MPM2002(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = -43.6977 - 0.0129037 * sal + 1.364e-4 * sal^2 +
    2885.378 / TK + 7.045159 * log(TK)

    K1 = 10.0^(-pK1) 

    pK2 = (-452.0940 + 13.142162 * sal - 8.101e-4 * sal^2 +
    21263.61 / TK + 68.483143 * log(TK) + (-581.4428 * sal +
    0.259601 * sal^2) / TK - 1.967035 * sal * log(TK))

    K2 = 10.0^(-pK2)

    return (; K1=K1, K2=K2)
end


# 12. Millero, et al., 2002: https://doi.org/10.1016/S0967-0637(02)00093-6
# Salinity between 34-37 PSU, temperature between -1.6-35 C
# pH on seawater scale, intended for field measurements
# NOTES FROM CO2SYS:
# Calculated from overdetermined WOCE-era field measurements
# sigma for pK1 is reported to be 0.005
# sigma for pK2 is reported to be 0.008
# This is from page 1715
function Millero2002(; temp_c, sal, kwargs...)

    pK1 = 6.359 - 0.00664 * sal - 0.01322 * temp_c + 4.989e-5 * temp_c^2

    K1 = 10.0^(-pK1) 

    pK2 = 9.867 - 0.01314 * sal - 0.01904 * temp_c + 2.448e-5 * temp_c^2

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end 


# 13. Millero, et al., 2006: https://doi.org/10.1016/j.marchem.2005.12.001
# Salinity between 1-50 PSU, temperature between 0-50 C
# pH on seawater scale, real seawater
function Millero2006(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    pK2₀ = -90.18333 + 5143.692 / TK + 14.613358 * log(TK)
    A₁ = 13.4191 * sqrt(sal) + 0.0331 * sal - 5.33e-5 * sal^2
    B₁ = -530.123 * sqrt(sal) - 6.103 * sal
    C₁ = -2.06950 * sqrt(sal)
    A₂ = 21.0894 * sqrt(sal) + 0.1248 * sal - 3.687e-4 * sal^2
    B₂ = -772.483 * sqrt(sal) - 20.051 * sal
    C₂ = -3.3336 * sqrt(sal)

    pK1 = A₁ + B₁ / TK + C₁ * log(TK) + pK1₀

    K1 = 10.0^(-pK1) 

    pK2 = A₂ + B₂ / TK + C₂ * log(TK) + pK2₀

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end 


# 14. Millero, et al., 2010: https://www.cabidigitallibrary.org/doi/full/10.5555/20103100380
# Salinity between 1-50 PSU, temperature between 0-50 C
# pH on seawater scale, real seawater
# For estuarine use
function Millero2010(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    pK2₀ = -90.18333 + 5143.692 / TK + 14.613358 * log(TK)
    A₁ = 13.4038 * sqrt(sal) + 0.03206 * sal - 5.242e-5 * sal^2
    B₁ = -530.659 * sqrt(sal) - 5.8210 * sal
    C₁ = -2.0664 * sqrt(sal)
    A₂ = 21.3728 * sqrt(sal) + 0.1218 * sal - 3.688e-4 * sal^2
    B₂ = -788.289 * sqrt(sal) - 19.189 * sal
    C₂ = -3.374 * sqrt(sal)

    pK1 = pK1₀ + A₁ + B₁ / TK + C₁ * log(TK)

    K1 = 10.0^(-pK1)

    pK2 = pK2₀ + A₂ + B₂ / TK + C₂ * log(TK)

    K2 = 10.0^(-pK2)

    return (; K1=K1, K2=K2)
end 


# 15. Waters, et al., 2014: https://doi.org/10.1016/j.marchem.2014.07.004
# Salinity between 1-50 PSU, temperature between 0-50 C
# pH on seawater scale, real seawater
# NOTES FROM CO2SYS:
# Corrigendum to "The free proton concentration scale for seawater pH".
# Effectively, this is an update of Millero (2010) formulation
# (WhichKs==14)
# Constants for K's on the SWS;
function Waters2014(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    pK2₀ = -90.18333 + 5143.692 / TK + 14.613358 * log(TK)
    A₁ = 13.409160 * sqrt(sal) + 0.031646 * sal - 5.1895e-5 * sal^2
    B₁ = -531.3642 * sqrt(sal) - 5.713 * sal
    C₁ = -2.0669166 * sqrt(sal)
    A₂ = 21.225890 * sqrt(sal) + 0.12450870 * sal - 3.7243e-4 * sal^2
    B₂ = -779.3444 * sqrt(sal) - 19.91739 * sal
    C₂ = -3.3534679 * sqrt(sal)

    pK1 = pK1₀ + A₁ + B₁ / TK + C₁ * log(TK)

    K1 = 10.0^(-pK1)

    pK2 = pK2₀ + A₂ + B₂ / TK + C₂ * log(TK)

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end


# 16. Schockman & Byrne, 2021: https://doi.org/10.1016/j.gca.2021.02.008
# Salinity between 19-41 PSU, temperature between 15-35 C (CO2SYS v3.2 option 17)
# pH on total scale, converted to SWS if /SWStoTOT is included
function SB2020(; temp_c, sal, ST=nothing, FT=nothing, KS=nothing, KF=nothing, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin
    
    # K1 Reference: Waters, Millero, and Woosley (2014) coefficients used by MATLAB
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    A₁   = 13.568513 * sqrt(sal) + 0.031645 * sal - 5.3834e-5 * sal^2
    B₁   = -539.2304 * sqrt(sal) - 5.635 * sal
    C₁   = -2.0901396 * sqrt(sal)

    pK1 = pK1₀ + A₁ + B₁ / TK + C₁ * log(TK)
    K1 = 10.0^(-pK1) 
    
    # K2 Reference: Schockman & Byrne (2021) — Checked: These perfectly match MATLAB!
    e1 = 116.8067
    e2 = -3655.02
    e3 = -16.45817
    e4 = 0.04523
    e5 = -0.615
    e6 = -0.0002799
    e7 = 4.969

    pK2 = e1 + e2 / TK + e3 * log(TK) + e4 * sal + e5 * sqrt(sal) +
          e6 * sal^2 + e7 * (sal / TK)

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end


# 17. Papadimitriou, et al., 2018: https://doi.org/10.1016/j.gca.2017.09.037
# Salinity between 33-100 PSU, temperature between (-6)-25 C (the paper's abstract)
# pH on total scale, converted to SWS if converted to SWS if /SWStoTOT is included
function Papadimitriou2018(; temp_c, sal, ST, FT, KS, KF, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = -176.48 + 6.14528 * sqrt.(sal) - 0.127714 * sal + 7.396 * 
    10.0^(-5) * sal^2 + (9914.37 - 622.886 * sqrt(sal) + 29.714 * sal) /
    TK + (26.05129 - 0.666812 * sqrt(sal)) * log(TK)

    pK2 = -323.52692 + 27.557655 * sqrt(sal) + 0.154922 * sal -
    2.48396 * 10.0^(-4) * sal^2 + (14763.287 - 1014.819 * sqrt(sal) -
    14.35223 * sal) / TK + (50.385807 - 4.4630415 * sqrt(sal)) * log(TK)

    K1 = 10.0^(-pK1) # / SWStoTOT(ST = ST, FT=FT, KS=KS, KF=KF) 

    K2 = 10.0^(-pK2) # / SWStoTOT(ST = ST, FT=FT, KS=KS, KF=KF) # converted to SWS 

    return (; K1=K1, K2=K2)
end


# 18. Sulpis, et al., 2020: https://doi.org/10.5194/os-16-847-2020
# Salinity between 31-38 PSU, temperature between -1.7-32 C (field data; CO2SYS v3.2 option 16)
# Total pH scale, real seawater (converted to SWS converted to SWS if /SWStoTOT is included)
function Sulpis2020(; temp_c, sal, ST, FT, KS, KF, kwargs...)
    TK = temp_c + 273.15     # Convert to Kelvin

    pK1 = 8510.63 / TK - 172.4493 + 26.32996 * log(TK) - 0.011555 * 
    sal + 0.0001152 * sal^2

    K1 = 10.0^(-pK1) # / SWStoTOT(ST=ST, FT=FT, KS=KS, KF=KF)

    pK2 =  4226.23 / TK - 59.4636 + 9.60817 * log(TK) - 0.01781 *
    sal + 0.0001122 * sal^2

    K2 = 10.0^(-pK2) # / SWStoTOT(ST=ST, FT=FT, KS=KS, KF=KF)

    return (; K1=K1, K2=K2)
end


# Calculations for total boron

# Case 8. Millero, 1979: https://doi.org/10.1016/0016-7037(79)90184-4
# Salinity of 0 PSU (freshwater), temperature between 0-50 C
function case8_BT(; kwargs...)
    return (; BT = 0.0) # in mol/kg-SW
end


# Case 6/7. Mehrbach 1973: https://doi.org/10.4319/lo.1973.18.6.0897
# Salinity between 19-43 PSU, temperature between 2-35 C
# NBS pH scale converted to SWS, real seawater
# NOTES FROM CO2SYS:
# This is .00001173*Sali, about 1% lower than Uppstrom's value
# Culkin, F., in Chemical Oceanography, ed. Riley and Skirrow, 1965:
# GEOSECS references this, but this value is not explicitly given here
function case67_BT(; sal, kwargs...)
    BT = 0.0004106 * sal / 35.0 # in mol/kg-SW

    # NOTES FROM CO2SYS:
    # this is .00001173*Sali
    # this is about 1# lower than Uppstrom's value
    # Culkin, F., in Chemical Oceanography,
    # ed. Riley and Skirrow, 1965:
    # GEOSECS references this, but this value is not explicitly given here
    return (; BT=BT) # in mol/kg-SW
end


# Calculating for other cases depends on which combination of KSO4 and BT desired

# Uppstrom, L., Deep-Sea Research 21:161-162, 1974
# NOTES FROM CO2SYS:
# Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
# this is .000416*Sali/35. = .0000119*Sali
# total_borate[FF] = (0.000232/10.811)*(Sal[FF]/1.80655); in mol/kg-SW.
function Uppstrom_BT(; sal, kwargs...)
    BT = 0.0004157 * sal / 35.0 # in mol/kg-SW
    return (; BT=BT) # in mol/kg-SW
end

# Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.
# Geochimica Et Cosmochimica Acta 74 (6): 1801â€“1811.
function Lee_BT(; sal, kwargs...)
    BT = 0.0004326 * sal / 35.0 # in mol/kg-SW

    return (; BT=BT) # in mol/kg-SW
end

# NOTES FROM CO2SYS:
# Note that the reference provided an equation for μmol/kg-sw
# This function divides it by a factor of 1e6 to convert to mol/kg-sw
function KSK18_BT(; sal, kwargs...)
    BT = (10.838 * sal + 13.821) / 1e6

    # NamedTuple, like the other *_BT functions - the call site in calculate_constants
    # does `.BT`, so returning a bare number made BT_method="KSK18" throw a FieldError.
    return (; BT=BT) # in mol/kg-SW
end

# Calculate FT
# Riley, J. P., Deep-Sea Research 12:219-220, 1965
function  calc_FT(; sal, kwargs...)
    FT = (0.000067 / 18.998) * (sal / 1.80655)
    return (; FT=FT) # in mol/kg-SW
end

# Calculate ST
# Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966
function calc_ST(; sal, kwargs...)
    ST = (0.14 / 96.062) * (sal / 1.80655)
    return (; ST=ST) # in mol/kg-SW
end


# Calculate K0
# Weiss, R. F., Marine Chemistry 2:203-215, 1974
function calc_K0(; temp_c, sal, kwargs...)
    TempK100 = (temp_c + 273.15) / 100
    lnK0 = (-60.2409 + 93.4517 / TempK100 + 23.3585 * log(TempK100) + 
    sal * (0.023517 - 0.023656 * TempK100 + 0.0047036 * TempK100^2))

    K0 = exp(lnK0) # in mol/kg-SW/atm
    return (; K0=K0)

end


# Calculate KS

# Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
# NOTES FROM CO2SYS:
# The goodness of fit is .021.
# It was given in mol/kg-H2O. I convert it to mol/kg-SW.
# TYPO on p. 121: the constant e9 should be e8.
# This is from eqs 22 and 23 on p. 123, and Table 4 on p 121

function Dickson_KS(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    # Ionic strength on the molal scale
    IonS = 19.924 * sal / (1000.0 - 1.005 * sal)

    lnKS = (-4276.1 / TK + 141.328 - 23.093 * log(TK)) + 
           (-13856.0 / TK + 324.57 - 47.986 * log(TK)) * sqrt(IonS) + 
           (35474.0 / TK - 771.54 + 114.723 * log(TK)) * IonS + 
           (-2698.0 / TK) * IonS^1.5 + (1776.0 / TK) * IonS^2

    # Convert from mol/kg-H2O to mol/kg-SW
    KS = exp(lnKS) * (1.0 - 0.001005109 * sal) 
    return (; KS=KS)
end

# Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
# NOTES FROM CO2SYS:
# KS was found by titrations with a hydrogen electrode
# of artificial seawater containing sulfate (but without F)
# at 3 salinities from 20 to 45 and artificial seawater NOT
# containing sulfate (nor F) at 16 salinities from 15 to 45,
# both at temperatures from 5 to 40 deg C.
# KS is on the Free pH scale (inherently so).
# It was given in mol/kg-H2O. I convert it to mol/kg-SW.
# He finds log(beta) which = my pKS;
# his beta is an association constant.
# The rms error is .0021 in pKS, or about .5# in KS.
# This is equation 20 on p. 33:
function Khoo_KS(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15

    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
    IonS = 19.924 * sal / (1000 - 1.005 * sal)

    pKS = 647.59 / TK - 6.3451 + 0.019085 * TK - 0.5208 * sqrt(IonS)
    KS = (10.0^(-pKS)) * (1 - 0.001005109 * sal) # mol/kg-SW

    return (; KS=KS)

end

function WM13_KS(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    logKS0 = (562.69486 - 102.5154 * log(TK) - 0.0001117033 * TK^2 + 0.2477538 * TK -
    13273.76 / TK)

    logKSK0 = ((4.24666 - 0.152671 * TK + 0.0267059 * TK * log(TK) - 0.000042128 *
    TK^2) * sal^0.5 + (0.2542181 - 0.00509534 * TK + 0.00071589 * TK * log(TK) * sal) +
    (-0.00291179 + 0.0000209968 * TK) * sal^1.5 -0.0000403724 * sal^2)

    kSO4 = (1 - 0.001005109 * sal) * 10^(logKS0 + logKSK0)

    # NamedTuple, like Dickson_KS and Khoo_KS. Returning a bare number made the three
    # disagree, which the call site papered over with an `isa Number` test on every call.
    return (; KS=kSO4) # in mol/kg-SW
end

# Calculate KF

# Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979
function Dickson_KF(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15

    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
    IonS = 19.924 * sal / (1000 - 1.005 * sal)
    lnKF = 1590.2 / TK - 12.641 + 1.525 * sqrt(IonS)
    KF = (exp(lnKF) * (1.0 - 0.001005109 * sal)) # in mol/kg-SW
    return (; KF=KF)
end

# Perez & Fraga, 1987: https://doi.org/10.1016/0304-4203(87)90036-3
# Salinity 10-40 PSU and Temperature of 9-33 C
# NOTE FrOM PYCO2SYS:
# Note that this is not currently used or an option in CO2SYS,
# despite the equations below appearing in CO2SYS.m (commented out).
function Perez_KF(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15

    lnKF = 874 / TK - 9.68 + 0.111 * sqrt(sal)
    KF = (exp(lnKF) * (1.0 - 0.001005109 * sal)) # # in mol/kg-SW
    return (; KF=KF)
end



# Calculate fH

# Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
"""
    calc_fH(temp_c, sal)

Activity coefficient ratio linking the NBS and seawater pH scales, `a(H)_NBS / m(H)_SWS`.

Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80), via CO2SYS.

Temperature is in Celsius, as everywhere else in the package; the polynomial itself is
fitted in Kelvin and the conversion happens here, so there is one place to get it wrong
rather than one per call site. Expect ~0.71 for seawater.
"""
function calc_fH(temp_c, sal)
    temp_k = temp_c + 273.15
    a, b, c, d = (1.2948, -2.036e-3, 4.607e-4, -1.475e-6)
    return a + b * temp_k + (c + d * temp_k) * sal^2
end


"`calc_fH` in the keyword form the `fH_method` registry calls."
Takahashi_fH(; temp_c, sal, kwargs...) = (; fH = calc_fH(temp_c, sal))

# Case #8: Millero, 1979: https://doi.org/10.1016/0016-7037(79)90184-4
# Salinity of 0 PSU (freshwater), temperature between 0-50 C
# Pure water, where the NBS and seawater scales coincide.
function case8_fH(; kwargs...)
    return (; fH = 1.0)
end

# Case #7: Mehrbach 1973: https://doi.org/10.4319/lo.1973.18.6.0897
# Salinity between 19-43 PSU, temperature between 2-35 C
# NBS pH scale converted to SWS, real seawater
function case7_fH(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    fH = (1.29 - 0.00204 * TK + (0.00046 - 0.00000148 * TK) * sal * sal)
    
    # NOTES FROM CO2SYS
    # Peng et al, Tellus 39B:439-458, 1987:
    # They reference the GEOSECS report, but round the value
    # given there off so that it is about .008 (1#) lower. It
    # doesn't agree with the check value they give on p. 456.
    return (; fH=fH)
end

# Calculate KB

# Dickson, A. G., Deep-Sea Research 37:755-766, 1990
function calc_KB(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    lnKBnum = (-8966.9 - 2890.53 * sqrt(sal) - 77.942 * sal + 1.728 *
    sqrt(sal) * sal - 0.0996 * sal^2)
    lnKB = (lnKBnum / TK + 148.0248 + 137.1942 * sqrt(sal) + 1.62142 *
    sal + (-24.4344 - 25.085 * sqrt(sal) - 0.2474 * sal) * log(TK) +
    0.053105 * sqrt(sal) * TK)

    KB = (exp(lnKB)) # native total scale
    return (; KB=KB)
end


# Calculate KW

# Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995
# At salinity 0 this is identical to Millero (1979)'s pure-water KW.
function calc_KW(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    lnKW = (148.9802 - 13847.26 / TK - 23.6521 * log(TK) + (-5.977 + 118.67 /
    TK + 1.0495 * log(TK)) * sqrt(sal) - 0.01615 * sal)
    KW = exp(lnKW) # native SWS scale
    return (; KW=KW)
end


# Calculate KP1, KP2, KP3, and KSi

# Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
# KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
# KSi was given on the SWS pH scale in molal units.
function calc_KP(; temp_c, sal, fH, kwargs...)
    TK = temp_c + 273.15
    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
    IonS = 19.924 * sal / (1000 - 1.005 * sal)

    lnKP1 = (-4576.752 / TK + 115.54 - 18.453 * log(TK) + (-106.736 / TK +
    0.69171) * sqrt(sal) + (-0.65643 / TK - 0.01844) * sal)
    KP1 = exp(lnKP1)

    lnKP2 = (-8814.715 / TK + 172.1033 - 27.927 * log(TK) + (-160.34 / TK +
    1.3566) * sqrt(sal) + (0.37335 / TK - 0.05778) * sal)
    KP2 = exp(lnKP2)

    lnKP3 = (-3070.75 / TK - 18.126 + (17.27039 / TK + 2.81197) * 
    sqrt(sal) + (-44.99486 / TK - 0.09984) * sal)
    KP3 = exp(lnKP3)

    lnKSi = (-8904.2 / TK + 117.4 - 19.334 * log(TK) + (-458.79 / TK +
    3.5913) * sqrt(IonS) + (188.74 / TK - 1.5998) * IonS + (-12.1652 /
    TK + 0.07871) * IonS^2)
    KSi = (exp(lnKSi) * (1 - 0.001005109 * sal))

    return (; KP1=KP1, KP2=KP2, KP3=KP3, KSi=KSi)
end



# Calculate KspA and KspC:
# Mucci, 1983: https://doi.org/10.2475/ajs.283.7.780
function calc_KspA(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15

    log10KspA = -171.945 - 0.077993 * TK + (2903.293 / TK) + 71.595 *
    log10(TK) + (-0.068393 + 0.0017276 * TK + (88.135 / TK)) * sqrt(sal) -
    0.10018 * sal + 0.0059415 * sal^(1.5)

    KspA = 10.0^(log10KspA)

    return (; KspA=KspA)
end 


function calc_KspC(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15

    log10KspC = -171.9065 - 0.077993 * TK + (2839.319 / TK) + 71.595 *
    log10(TK) + (-0.77712 + 0.0028426 * TK + (178.34 / TK)) * sqrt(sal) -
    0.07711 * sal + 0.0041249 * sal^(1.5)

    KspC = 10.0^(log10KspC)

    return (; KspC=KspC)
end


# Calculate KH2S

# NOTES FROM CO2SYS:
 # H2S  Millero et. al.( 1988)  Limnol. Oceanogr. 33,269-274.
# Yao and Millero, Aquatic Geochemistry 1:53-88, 1995. Total Scale.
# Yao Millero say equations have been refitted to SWS scale but not true as
# they agree with Millero 1988 which are on Total Scale.
# Also, calculations agree at high H2S with AquaEnv when assuming it is on
# Total Scale.
function calc_KH2S(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    lnKH2S = 225.838 - 13275.3 / TK - 34.6435 * log(TK) + 0.3449 * sqrt(sal) - 0.0247 * sal
    KH2S = exp(lnKH2S)
    return (; KH2S=KH2S)
end


# Calculate KNH3

# NOTES FROM CO2SYS:
# Yao and Millero, Aquatic Geochemistry 1:53-88, 1995   SWS
function Millero_KNH3(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    lnKNH3 = -6285.33 / TK + 0.0001635 * TK - 0.25444 + (0.46532 - 123.7184 / TK) *
    sqrt(sal) + (-0.01992 + 3.17556 / TK) * sal
    KNH3 = exp(lnKNH3)
    return (; KNH3=KNH3)
end


# NOTES FROM CO2SYS:
# Clegg Whitfield 1995
# Geochimica et Cosmochimica Acta, Vol. 59, No. 12. pp. 2403-2421
# eq (18)  Total scale   t=[-2 to 40 oC]  S=[0 to 40 ppt]   pK=+-0.00015
function Clegg_KNH3(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    pK = 9.244605 - 2729.33 * (1 / 298.15 - 1 / TK)
    pK += (0.04203362 - 11.24742 / TK) * sal^0.25
    pK += (
        -13.6416 + 1.176949 * sqrt(TK) - 0.02860785 * TK + 545.4834 / TK
    ) * sqrt(sal)
    pK += (
        -0.1462507 + 0.0090226468 * sqrt(TK) - 0.0001471361 * TK+ 10.5425 / TK
    ) * sal^1.5
    pK += (
        0.004669309 - 0.0001691742 * sqrt(TK) - 0.5677934 / TK
    ) * sal^2
    pK += (-2.354039e-05 + 0.009698623 / TK) * sal^2.5
    KNH3 = 10.0^-pK
    # 3. Convert from mol/kg-H2O to mol/kg-SW
    KNH3 = KNH3 * (1 - 0.001005109 * sal)

    return (; KNH3=KNH3)
end


# Calculate Ca

# NOTES FROM CO2SYS
# *** CalculateCaforGEOSECS:
# Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
# (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
# Culkin gives Ca = (.0213/40.078)*(Sal/1.80655) in mol/kg-SW
# which corresponds to Ca = .01030*Sal/35.
function Culkin_Ca(; sal, kwargs...)
    return (0.01026 * sal / 35) # in mol/kg-SW
end


function RT67_Ca(; sal, kwargs...)
    return (0.02128 / 40.087 * sal / 1.80655) # in mol/kg-SW
end


# --- Pressure corrections ------------------------------------------------------------------

"""
Millero 1995 volume and compressibility coefficients, one row per constant.

Each row holds the two polynomials in temperature that [`_pressure_factor`](@ref) evaluates:
`ΔV = a + b·T + c·T²` and `κ = (d + e·T + f·T²)/1000`. Data rather than code, because the
correction is the same three lines of arithmetic for all sixteen constants and only these
numbers differ.

Every row uses the same ÷1000 convention, KspA and KspC included. Folding that factor into
their literals instead (`-11.76e-3 + 0.3692e-3·T`) rounds fractionally differently and shifts
their pressure correction by an ulp, so keep the convention uniform.
"""
const PRESSURE_COEFFICIENTS = (
    K1   = (ΔV = (-25.5,  0.1271,  0.0),       κ = (-3.08,  0.0877, 0.0)),
    K2   = (ΔV = (-15.82, -0.0219, 0.0),       κ = (1.13,   -0.1475, 0.0)),
    KB   = (ΔV = (-29.48, 0.1622,  -0.002608), κ = (-2.84,  0.0,    0.0)),
    KW   = (ΔV = (-20.02, 0.1119,  -0.001409), κ = (-5.13,  0.0794, 0.0)),
    KF   = (ΔV = (-9.78,  -0.009,  -0.000942), κ = (-3.91,  0.054,  0.0)),
    KS   = (ΔV = (-18.03, 0.0466,  0.000316),  κ = (-4.53,  0.09,   0.0)),
    KP1  = (ΔV = (-14.51, 0.1211,  -0.000321), κ = (-2.67,  0.0427, 0.0)),
    KP2  = (ΔV = (-23.12, 0.1758,  -0.002647), κ = (-5.15,  0.09,   0.0)),
    KP3  = (ΔV = (-26.57, 0.202,   -0.003042), κ = (-4.08,  0.0714, 0.0)),
    KSi  = (ΔV = (-29.48, 0.1622,  -0.002608), κ = (-2.84,  0.0,    0.0)),
    KH2S = (ΔV = (-11.07, -0.009,  -0.000942), κ = (-2.89,  0.054,  0.0)),
    KNH3 = (ΔV = (-26.43, 0.0889,  -0.000905), κ = (-5.03,  0.0814, 0.0)),
    KspA = (ΔV = (-46.0,  0.5304,  0.0),       κ = (-11.76, 0.3692, 0.0)),
    KspC = (ΔV = (-48.76, 0.5304,  0.0),       κ = (-11.76, 0.3692, 0.0)),
)

"""
Millero 1983's freshwater coefficients, selected with `pressure_method = "Millero 1983"` and
published beside the `Millero 1979` K1/K2 parameterisation.

A table rather than a branch, so the special case is a few rows that `merge` overlays on
`PRESSURE_COEFFICIENTS`. KB is zeroed rather than omitted, which makes its correction a factor
of exactly 1 — immaterial either way, since BT is ~zero in freshwater.
"""
const MILLERO1979_PRESSURE = (
    K1 = (ΔV = (-30.54, 0.1849, -0.0023366), κ = (-5.74, 0.093,  -0.001896)),
    K2 = (ΔV = (-29.81, 0.115,  -0.001816),  κ = (-5.74, 0.093,  -0.001896)),
    KW = (ΔV = (-25.6,  0.2324, -0.0036246), κ = (-7.33, 0.1368, -0.001233)),
    KB = (ΔV = (0.0,    0.0,    0.0),        κ = (0.0,   0.0,    0.0)),
)

"""
    _pressure_factor(temp_c, pres_bar, coefficients)

What an equilibrium constant is multiplied by to move it from the surface to `pres_bar`.
"""
function _pressure_factor(temp_c, pres_bar, coefficients)
    ΔV = coefficients.ΔV[1] + coefficients.ΔV[2] * temp_c + coefficients.ΔV[3] * temp_c^2
    κ = (coefficients.κ[1] + coefficients.κ[2] * temp_c + coefficients.κ[3] * temp_c^2) / 1000
    RT = GAS_CONSTANT * (temp_c + 273.15)
    return exp((-ΔV + 0.5 * κ * pres_bar) * pres_bar / RT)
end


# --- Option lookup tables -----------------------------------------------------------------

"""
Modern seawater calcium, as a parameterisation so that `Ca_method` has interchangeable options
rather than several plus a fallback. Independent of salinity by definition: this is the
composition Kgen assumes for the KGen correction, not a concentration to be diluted.
"""
Modern_Ca(; kwargs...) = MODERN_CALCIUM

"The standard pressure coefficients, in the keyword form the `pressure_method` registry calls."
Millero1995_pressure(; kwargs...) = PRESSURE_COEFFICIENTS

"""
Millero (1983)'s freshwater coefficients for K1, K2 and KW over the standard table.

Identical field names and types to [`Millero1995_pressure`](@ref)'s, so the choice stays
inferable — a `Union` of two differently-shaped tables would leave every corrected constant
untyped.
"""
Millero1983_pressure(; kwargs...) = merge(PRESSURE_COEFFICIENTS, MILLERO1979_PRESSURE)

"""
How much of each dissolved species the water holds, by the name a [`KSettings`](@ref) field
selects.

Composition, not speciation: these parameterisations give a *concentration* from salinity.

Rows are `(fn = …)`, plus `intended_for` on an option that only makes sense beside a
particular K1/K2 parameterisation, which [`KSettings`](@ref) warns about when it is used with
any other. See [`SPECIATION`](@ref) for what the shared `default`/`options` shape is for.
"""
const COMPOSITION = (
    BT_method = (
        default = "Uppstrom",
        options = Dict(
            "Uppstrom"   => (fn = Uppstrom_BT,),
            "Lee"        => (fn = Lee_BT,),
            "KSK18"      => (fn = KSK18_BT,),
            # GEOSECS and Peng et al.'s boron, 0.0004106·S/35, about 1 % below Uppstrom.
            "Culkin"     => (fn = case67_BT,),
            # No boron at all: the freshwater set that goes with Millero 1979.
            "freshwater" => (fn = case8_BT, intended_for = ("Millero 1979",))),
        ),
    Ca_method = (
        default = "modern",
        options = Dict(
            "modern" => (fn = Modern_Ca,),
            "Culkin" => (fn = Culkin_Ca,),
            "RT67"   => (fn = RT67_Ca,))
        ),
)

"""
The equilibrium constants, by the name a [`KSettings`](@ref) field selects, each with the pH
scale it is fitted on.

Rows are `(fn = …, scale = …)`. `K1_K2_method` rows also carry:

- `temp_c`, `sal`: the `(low, high)` range the parameterisation was fitted over, in °C and
  practical salinity. A calculation outside it warns; see `_check_K_settings`.
- `source`: where that range comes from, for the warning and the README.
- `companions`: the sets of other constants the parameterisation was published or used with,
  as `(name, settings)` pairs. [`KSettings`](@ref) warns when none of them is selected. Empty
  for the parameterisations that are self-contained.

`:free` on KS and KF is recorded but never acted on: the free-scale values are exactly what
`SWStoTOT` and `FREEtoTOT` are built from. It is here so the assumption sits with the
parameterisation rather than in a comment beside a pressure correction.

The call site for each argument is shared, so every row within one entry must take the same
arguments and return the same shape. That works because every parameterisation in this file
ends in `kwargs...`: one argument set calls all of them and each takes what it needs. All of
them return a NamedTuple.
"""
const SPECIATION = (
    KSO4_method = (
        default = "Dickson",
        options = Dict(
            "Dickson" => (fn = Dickson_KS, scale = :free),
            "Khoo" => (fn = Khoo_KS, scale = :free),
            "WM13" => (fn = WM13_KS, scale = :free))
        ),
    KF_method   = (
        default = "Dickson",
        options = Dict(
            "Dickson" => (fn = Dickson_KF, scale = :free),
            "Perez" => (fn = Perez_KF, scale = :free))
        ),
    KNH3_method = (
        default = "Millero",
        options = Dict(
            "Millero" => (fn = Millero_KNH3, scale = :sws),
            "Clegg" => (fn = Clegg_KNH3, scale = :total))
        ),
    # Ranges are from the CO2SYS v3.2 header unless `source` says otherwise. CO2SYS option
    # numbers are given where a range came from there, so it can be checked against the same
    # table.
    K1_K2_method = (
        default = "KGen",
        options = Dict(
            "Roy 1993" => (fn = Roy1993, scale = :total,
                temp_c = (0.0, 45.0), sal = (5.0, 45.0),
                source = "Roy et al. (1993); CO2SYS option 1", companions = ()),
            "GP 1989" => (fn = GP1989, scale = :sws,
                temp_c = (-1.0, 40.0), sal = (10.0, 50.0),
                source = "Goyet & Poisson (1989); CO2SYS option 2", companions = ()),
            "Hansson 1973" => (fn = Hansson1973, scale = :sws,
                temp_c = (2.0, 35.0), sal = (20.0, 40.0),
                source = "Hansson (1973) refit by Dickson & Millero (1987); CO2SYS option 3",
                companions = ()),
            "DM 1987" => (fn = DM1987, scale = :sws,
                temp_c = (2.0, 35.0), sal = (20.0, 40.0),
                source = "Mehrbach et al. (1973) refit by Dickson & Millero (1987); " *
                         "CO2SYS option 4", companions = ()),
            "HM 1973" => (fn = HM1973, scale = :sws,
                temp_c = (2.0, 35.0), sal = (20.0, 40.0),
                source = "Hansson and Mehrbach combined by Dickson & Millero (1987); " *
                         "CO2SYS option 5", companions = ()),
            # One function, two historical sets of companions. Each used to be its own
            # K_method ("Mehrbach 1973 A"/"B"), which switched the companions silently.
            "Mehrbach 1973" => (fn = Mehrbach1973, scale = :sws,
                temp_c = (2.0, 35.0), sal = (19.0, 43.0),
                source = "Mehrbach et al. (1973), as used by GEOSECS and Peng et al.; " *
                         "CO2SYS options 6 and 7",
                companions = (
                    (name = "GEOSECS (Takahashi et al., 1982; CO2SYS option 6)",
                     settings = (fH_method = "Takahashi 1982", BT_method = "Culkin",
                                 Ca_method = "Culkin")),
                    (name = "Peng et al. (1987; CO2SYS option 7)",
                     settings = (fH_method = "Peng 1987", BT_method = "Culkin",
                                 Ca_method = "Culkin")),
                )),
            "Millero 1979" => (fn = Millero1979, scale = :nbs,
                temp_c = (0.0, 50.0), sal = (0.0, 0.0),
                source = "Millero (1979), pure water only; CO2SYS option 8",
                companions = (
                    (name = "freshwater (CO2SYS option 8)",
                     settings = (fH_method = "freshwater", BT_method = "freshwater",
                                 pressure_method = "Millero 1983")),
                )),
            "CW 2003" => (fn = CW2003, scale = :sws,
                temp_c = (2.0, 35.0), sal = (0.0, 49.0),
                source = "Cai & Wang (1998); CO2SYS option 9", companions = ()),
            "Lueker 2000" => (fn = Lueker2000, scale = :total,
                temp_c = (2.0, 35.0), sal = (19.0, 43.0),
                source = "Lueker et al. (2000); CO2SYS option 10", companions = ()),
            "MPM 2002" => (fn = MPM2002, scale = :sws,
                temp_c = (0.0, 45.0), sal = (5.0, 42.0),
                source = "Mojica Prieto & Millero (2002); CO2SYS option 11",
                companions = ()),
            "Millero 2002" => (fn = Millero2002, scale = :sws,
                temp_c = (-1.6, 35.0), sal = (34.0, 37.0),
                source = "Millero et al. (2002); CO2SYS option 12", companions = ()),
            "Millero 2006" => (fn = Millero2006, scale = :sws,
                temp_c = (0.0, 50.0), sal = (1.0, 50.0),
                source = "Millero et al. (2006); CO2SYS option 13", companions = ()),
            "Millero 2010" => (fn = Millero2010, scale = :sws,
                temp_c = (0.0, 50.0), sal = (1.0, 50.0),
                source = "Millero (2010); CO2SYS option 14", companions = ()),
            "Waters 2014" => (fn = Waters2014, scale = :sws,
                temp_c = (0.0, 50.0), sal = (1.0, 50.0),
                source = "Waters et al. (2014); CO2SYS option 15", companions = ()),
            "Sulpis 2020" => (fn = Sulpis2020, scale = :total,
                temp_c = (-1.7, 32.0), sal = (31.0, 38.0),
                source = "Sulpis et al. (2020); CO2SYS option 16", companions = ()),
            "SB 2020" => (fn = SB2020, scale = :total,
                temp_c = (15.0, 35.0), sal = (19.0, 41.0),
                source = "Schockman & Byrne (2021); CO2SYS option 17", companions = ()),
            "Papadimitriou 2018" => (fn = Papadimitriou2018, scale = :total,
                temp_c = (-6.0, 25.0), sal = (33.0, 100.0),
                source = "Papadimitriou et al. (2018), the paper's abstract",
                companions = ()),
            # No K1/K2 function: Kgen hands back a whole bundle, already pressure-corrected
            # and on the total scale, assembled separately in `calculate_constants`.
            #
            # The range is the one this package's README has always given for KGen. Kgen's
            # own documentation states none, and its modern-seawater K1/K2 follow Dickson,
            # Sabine & Christian (2007), i.e. Lueker et al. (2000), fitted over 2-35 °C and
            # S 19-43 — so this is worth confirming with Kgen's authors.
            "KGen" => (fn = nothing, scale = :total,
                temp_c = (-2.0, 50.0), sal = (0.0, 50.0),
                source = "Kgen (Whiteford et al., 2025), range as given in the README",
                companions = ()),
        )
    ),
)

"""
How measured pH and pressure relate to the constants, by the name a [`KSettings`](@ref) field
selects.

`fH_method` is the activity coefficient of H⁺ that relates the NBS pH scale to the seawater
scale. It converts NBS-fitted constants (Mehrbach 1973, Millero 1979) and `pHNBS`.

`pressure_method` is the table of partial molal volumes and compressibilities the constants
are corrected to pressure with. Under KGen it reaches only KH2S and KNH3, because Kgen applies
its own correction to the rest.
"""
const CONVERSION = (
    fH_method = (
        default = "Takahashi 1982",
        options = Dict(
            "Takahashi 1982" => (fn = Takahashi_fH,),
            "Peng 1987"      => (fn = case7_fH, intended_for = ("Mehrbach 1973",)),
            "freshwater"     => (fn = case8_fH, intended_for = ("Millero 1979",))),
        ),
    pressure_method = (
        default = "Millero 1995",
        options = Dict(
            "Millero 1995" => (fn = Millero1995_pressure,),
            "Millero 1983" => (fn = Millero1983_pressure, intended_for = ("Millero 1979",))),
        ),
)

"""
How Kgen corrects the constants for seawater composition. Not a parameterisation, so it has
no function, but it is chosen in the same [`KSettings`](@ref) and validated the same way.

Only the polynomial approximation is implemented in Kgen.jl; the full MyAMI model is Python
only, so `"calculate"` is refused by name rather than listed.
"""
const KGEN_OPTIONS = (
    MyAMI_mode = (default = "approximate", options = Dict("approximate" => (fn = nothing,))),
)

"""
Raise for a method name that names no parameterisation.

`valid` is built from the table when the dispatch below is generated, so the list of names
cannot drift from the ones actually accepted, and the error path touches no table at all.
"""
_unknown_method(argument::Symbol, name, valid::String) =
    throw(ArgumentError("Unknown $argument: \"$name\".\nChoose one of: $valid"))

"The 'you could have meant' list for one argument, built at generation time."
_valid_names(entry) =
    join(sort!(collect(keys(entry.options))), ", ") * ". \"default\" is $(entry.default)."

"""
    _apply(Val(:KSO4_method), name; kwargs...)

Call the parameterisation `name` selects for one setting, resolving `"default"` and rejecting
anything the table does not list.

**Generated from [`COMPOSITION`](@ref), [`SPECIATION`](@ref) and [`CONVERSION`](@ref), not a
lookup into them.** A `Dict` holding several different functions can only be typed
`Dict{String,Function}`, so calling through it is a dynamic dispatch returning `Any` — and
that `Any` spreads into every constant derived from the result, leaving the whole bundle
untyped. Interpolating the function objects into a comparison chain at load time keeps every
call statically known while the tables stay the single place a parameterisation is
registered. Do not replace this with a `Dict` lookup - it is much slower.
"""
function _apply end

"""
    _native_scale(Val(:K1_K2_method), name) -> :total | :sws | :nbs

The pH scale the parameterisation `name` selects is fitted on.

Generated from [`SPECIATION`](@ref) only — the composition parameterisations give
concentrations, which are not on a pH scale.

Resolves `"default"` exactly as [`_apply`](@ref) does, which is the point of going through it
at all: a scale test written against the raw argument, such as
`KNH3_method == "Millero" ? sws_factor : 1.0`, silently misses the `"default"` that selects
Millero, and leaves a seawater-scale constant in a total-scale bundle.
"""
function _native_scale end

for table in (COMPOSITION, SPECIATION, CONVERSION), (argument, entry) in pairs(table)
    body = Expr(:block)
    push!(body.args, :(name == "default" && (name = $(entry.default))))

    for (option, row) in entry.options
        isnothing(row.fn) && continue    # KGen, which has no K1/K2 function
        push!(body.args, :(name == $option && return $(row.fn)(; kwargs...)))
    end

    @eval function _apply(::Val{$(QuoteNode(argument))}, name; kwargs...)
        $body
        _unknown_method($(QuoteNode(argument)), name, $(_valid_names(entry)))
    end
end

for (argument, entry) in pairs(SPECIATION)
    body = Expr(:block)
    push!(body.args, :(name == "default" && (name = $(entry.default))))

    for (option, row) in entry.options
        push!(body.args, :(name == $option && return $(QuoteNode(row.scale))))
    end

    @eval function _native_scale(::Val{$(QuoteNode(argument))}, name)
        $body
        _unknown_method($(QuoteNode(argument)), name, $(_valid_names(entry)))
    end
end

"""
The fitted range of every K1/K2 parameterisation, keyed by name, for the per-call check.

Built from [`SPECIATION`](@ref) at load, so the registry stays the one place a range is
written. A separate concretely-typed `Dict` because the registry's own rows differ in shape and
so only fit a `Dict{String, Any}`; looking a range up in that on every calculation would be a
dynamic dispatch on the hot path.
"""
const K1_K2_RANGES = Dict{String, @NamedTuple{temp_c::NTuple{2, Float64},
                                             sal::NTuple{2, Float64}, source::String}}(
    name => (temp_c = row.temp_c, sal = row.sal, source = row.source)
    for (name, row) in SPECIATION.K1_K2_method.options)


# --- Choosing the constants ---------------------------------------------------------------

"""
    KSettings(; K1_K2_method, KSO4_method, KF_method, KNH3_method, BT_method, Ca_method,
                fH_method, pressure_method, MyAMI_mode)

Which parameterisation to use for each equilibrium constant, total and conversion.

Name only the choices that differ from the defaults; everything else keeps its default.
Anywhere a `KSettings` is accepted, a NamedTuple of the same fields is too, so these are
equivalent:

```jldoctest
julia> carbon_system(TA = 2300.0, DIC = 2000.0,
                     K_settings = (K1_K2_method = "Roy 1993", KF_method = "Perez")).pHtot ==
       carbon_system(TA = 2300.0, DIC = 2000.0,
                     K_settings = KSettings(K1_K2_method = "Roy 1993",
                                            KF_method = "Perez")).pHtot
true
```

| field | options | default |
|---|---|---|
| `K1_K2_method` | `"KGen"`, `"Roy 1993"`, `"GP 1989"`, `"Hansson 1973"`, `"DM 1987"`, `"HM 1973"`, `"Mehrbach 1973"`, `"Millero 1979"`, `"CW 2003"`, `"Lueker 2000"`, `"MPM 2002"`, `"Millero 2002"`, `"Millero 2006"`, `"Millero 2010"`, `"Waters 2014"`, `"Sulpis 2020"`, `"SB 2020"`, `"Papadimitriou 2018"` | `"KGen"` |
| `KSO4_method` | `"Dickson"`, `"Khoo"`, `"WM13"` | `"Dickson"` |
| `KF_method` | `"Dickson"`, `"Perez"` | `"Dickson"` |
| `KNH3_method` | `"Millero"`, `"Clegg"` | `"Millero"` |
| `BT_method` | `"Uppstrom"`, `"Lee"`, `"KSK18"`, `"Culkin"`, `"freshwater"` | `"Uppstrom"` |
| `Ca_method` | `"modern"`, `"Culkin"`, `"RT67"` | `"modern"` |
| `fH_method` | `"Takahashi 1982"`, `"Peng 1987"`, `"freshwater"` | `"Takahashi 1982"` |
| `pressure_method` | `"Millero 1995"`, `"Millero 1983"` | `"Millero 1995"` |
| `MyAMI_mode` | `"approximate"` | `"approximate"` |

`"default"` names the default for any field, and a `Symbol` is read as its string.

Each choice is independent: nothing here changes another field. Where a K1/K2
parameterisation was published with particular companions — Mehrbach 1973 with the GEOSECS
or Peng et al. sets, Millero 1979 with the freshwater set — selecting it without them warns,
naming the sets and what is selected instead, and the calculation uses what you selected. The
same goes for a companion-only option (`fH_method = "Peng 1987"`, the freshwater options)
selected beside a different K1/K2 parameterisation, and for `KSO4_method`, `KF_method` or
`pressure_method` chosen alongside `"KGen"`, which computes those constants itself.

Whether the chosen K1/K2 parameterisation suits the temperature and salinity of the water is
checked per calculation instead, since that depends on the inputs; see the README for every
parameterisation's fitted range.

A supplied `BT`, `Ca` or `Mg` value overrides `BT_method` or `Ca_method`.

`KSettings(settings; changes...)` copies `settings` with some fields changed.
"""
struct KSettings
    K1_K2_method::String
    KSO4_method::String
    KF_method::String
    KNH3_method::String
    BT_method::String
    Ca_method::String
    fH_method::String
    pressure_method::String
    MyAMI_mode::String
end

"The registry entry a `KSettings` field is validated against."
@inline _setting_entry(field::Symbol) =
    haskey(SPECIATION, field)   ? SPECIATION[field] :
    haskey(COMPOSITION, field)  ? COMPOSITION[field] :
    haskey(CONVERSION, field)   ? CONVERSION[field] :
                                  KGEN_OPTIONS[field]

"What each `KSettings` field resolves to when not named."
const K_SETTINGS_DEFAULTS = NamedTuple{fieldnames(KSettings)}(
    map(f -> _setting_entry(f).default, fieldnames(KSettings)))

"""
Names that used to be accepted for a field, mapped to what to write now.

Read only to build an error, so the old names cost nothing until someone uses one. A retired
name is an error rather than a silent translation inside `KSettings`, because each of these
used to change *other* constants as well, and the point of `KSettings` is that nothing does.
"""
const RETIRED_K1_K2_NAMES = Dict(
    "Mehrbach 1973 A" =>
        "K1_K2_method = \"Mehrbach 1973\". For the GEOSECS set it used to imply, also set " *
        "BT_method = \"Culkin\" and Ca_method = \"Culkin\"",
    "Mehrbach 1973 B" =>
        "K1_K2_method = \"Mehrbach 1973\". For the Peng et al. set it used to imply, also " *
        "set fH_method = \"Peng 1987\", BT_method = \"Culkin\" and Ca_method = \"Culkin\"",
    "automatic" =>
        "a named K1_K2_method. \"automatic\" chose one from temperature and salinity, so the " *
        "constants changed between samples without saying so; each calculation now warns " *
        "instead when the chosen parameterisation does not suit the water. The README " *
        "tabulates every parameterisation's range",
)

"Resolve and validate one field's value."
function _resolve_setting(field::Symbol, value)
    value isa Union{AbstractString, Symbol} || throw(ArgumentError(
        "K_settings.$field must be a string naming a method, got $(repr(value))"))
    name = String(value)
    entry = _setting_entry(field)
    name == "default" && return entry.default
    haskey(entry.options, name) && return name

    field === :K1_K2_method && haskey(RETIRED_K1_K2_NAMES, name) && throw(ArgumentError(
        "K1_K2_method \"$name\" is no longer an option; use $(RETIRED_K1_K2_NAMES[name])."))
    field === :MyAMI_mode && name == "calculate" && throw(ArgumentError(
        "MyAMI_mode \"calculate\" is not available: Kgen.jl implements only the polynomial " *
        "approximation to MyAMI. The full model is in the Python Kgen."))
    _unknown_method(field, name, _valid_names(entry))
end

"""
Reject field names `KSettings` does not have, naming the likely intended one.

`K_method` gets its own message because it is the old name for the most common choice.
"""
function _reject_unknown_settings(names)
    fields = fieldnames(KSettings)
    unknown = [n for n in names if !(n in fields)]
    isempty(unknown) && return nothing

    squash(n) = lowercase(replace(String(n), "_" => ""))
    lines = map(unknown) do name
        name === :K_method && return "K_method is now K1_K2_method"
        near = [f for f in fields if squash(f) == squash(name)]
        isempty(near) ? "$name is not a K_settings field" : "$name (did you mean $(only(near))?)"
    end
    throw(ArgumentError("unrecognised K_settings field(s):\n  " * join(lines, "\n  ") *
                        "\nThe fields are: " * join(fields, ", ")))
end

function KSettings(; kwargs...)
    _reject_unknown_settings(keys(kwargs))
    values = map(f -> _resolve_setting(f, get(kwargs, f, "default")), fieldnames(KSettings))
    settings = KSettings(values...)
    _warn_companions(settings)
    return settings
end

KSettings(settings::KSettings; changes...) =
    KSettings(; merge(_non_default_fields(settings), NamedTuple(changes))...)

"""
The fields of `settings` that differ from the defaults, as a NamedTuple.

What `show` prints and what an explicitly-chosen option means for the warnings: a field set to
its default is indistinguishable from one left out, and should be.
"""
function _non_default_fields(settings::KSettings)
    names = Tuple(f for f in fieldnames(KSettings)
                  if getfield(settings, f) != K_SETTINGS_DEFAULTS[f])
    return NamedTuple{names}(map(f -> getfield(settings, f), names))
end

Base.:(==)(a::KSettings, b::KSettings) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(KSettings))
Base.hash(s::KSettings, h::UInt) =
    foldr((f, acc) -> hash(getfield(s, f), acc), fieldnames(KSettings); init = hash(:KSettings, h))

"Printed as the call that rebuilds it, naming only what differs from the defaults."
function Base.show(io::IO, settings::KSettings)
    changed = _non_default_fields(settings)
    print(io, "KSettings(", join(("$k = $(repr(v))" for (k, v) in pairs(changed)), ", "), ")")
end

"""
    _as_K_settings(x) -> KSettings

Accept a `KSettings`, a NamedTuple of its fields, or `nothing` for the defaults.
"""
_as_K_settings(settings::KSettings) = settings
_as_K_settings(::Nothing) = KSettings()
_as_K_settings(settings::NamedTuple) = KSettings(; settings...)
_as_K_settings(other) = throw(ArgumentError(
    "K_settings must be a NamedTuple such as (K1_K2_method = \"Roy 1993\",) or a " *
    "KSettings, got $(repr(other))"))

"`settings` as the literal a user would write, for warning messages."
_settings_literal(changes::NamedTuple) =
    "(" * join(("$k = $(repr(v))" for (k, v) in pairs(changes)), ", ") *
    (length(changes) == 1 ? ",)" : ")")

"""
Warn about choices that do not belong together. Runs when a `KSettings` is built, because
these depend only on the settings, not on the water.

Nothing is changed: the warnings say what the published method would have used, and the
calculation uses what was selected.

`maxlog = 1` with an `_id` per combination, so a solver built once and broadcast warns once,
and so does a keyword call repeated over a million rows.
"""
function _warn_companions(s::KSettings)
    method = s.K1_K2_method
    row = SPECIATION.K1_K2_method.options[method]

    # A K1/K2 parameterisation selected without any of its published companion sets.
    if !isempty(row.companions)
        if !any(set -> all(getfield(s, k) == v for (k, v) in pairs(set.settings)),
                row.companions)
            fields = unique(Iterators.flatten(keys(set.settings) for set in row.companions))
            sets = join(("  $(set.name): K_settings = " *
                         _settings_literal(merge((K1_K2_method = method,), set.settings))
                         for set in row.companions), "\n")
            current = join(("\"$(getfield(s, f))\" selected for $f" for f in fields), ", ")
            @warn "You've chosen \"$method\" for K1_K2_method. This K1/K2 method is " *
                  "intended to be used with " *
                  (length(row.companions) == 1 ? "this set" : "one of these sets") *
                  " of companion constants:\n$sets\nYou currently have $current. Your " *
                  "selections are used as given." maxlog=1 _id=(:companions, s)
        end
    end

    # An option that only makes sense beside one K1/K2 parameterisation, used with another.
    for field in (:BT_method, :fH_method, :pressure_method)
        option = getfield(s, field)
        intended = get(_setting_entry(field).options[option], :intended_for, nothing)
        if !isnothing(intended) && !(method in intended)
            @warn "You've chosen \"$option\" for $field. This option is intended to be " *
                  "used with K1_K2_method = " * join(("\"$m\"" for m in intended), " or ") *
                  ". You currently have \"$method\" selected for K1_K2_method." maxlog=1 _id=(
                      :intended_for, field, option, method)
        end
    end

    # Choices Kgen does not consult. Kgen computes KS, KF and the pressure correction of
    # everything in its bundle itself.
    if method == "KGen"
        for field in (:KSO4_method, :KF_method, :pressure_method)
            option = getfield(s, field)
            option == K_SETTINGS_DEFAULTS[field] && continue
            reach = field === :pressure_method ?
                "it only reaches KH2S and KNH3, because Kgen pressure-corrects the rest " *
                "of its constants itself" :
                "Kgen computes $(field === :KSO4_method ? "KSO4" : "KF") itself"
            @warn "You've chosen \"$option\" for $field, but you currently have \"KGen\" " *
                  "selected for K1_K2_method, so $reach. Choose a different K1_K2_method " *
                  "for $field to apply." maxlog=1 _id=(:kgen_ignores, field, option)
        end
    end
    return nothing
end


"The defaults, built once. `KSettings()` validates and checks companions on every call."
const DEFAULT_K_SETTINGS = KSettings()


# --- Settings given the old way -----------------------------------------------------------

"""
The keyword arguments that chose the constants before [`KSettings`](@ref), and which still
work, with a warning.
"""
const LEGACY_METHOD_KEYWORDS = (:K_method, :KSO4_method, :BT_method, :KF_method,
                                :KNH3_method, :Ca_method, :MyAMI_mode)

"""
What an old `K_method` value did, as `KSettings` fields.

The three that used to switch companions silently translate to those companions explicitly,
so an old script gets the numbers it always did.
"""
function _translate_K_method(name)
    name = String(name)
    name in ("default", "KGen") && return (K1_K2_method = "KGen",)
    name == "Mehrbach 1973 A" && return (K1_K2_method = "Mehrbach 1973", BT_method = "Culkin")
    name == "Mehrbach 1973 B" && return (K1_K2_method = "Mehrbach 1973",
                                         fH_method = "Peng 1987")
    name == "Millero 1979" && return (K1_K2_method = "Millero 1979",
                                      fH_method = "freshwater", BT_method = "freshwater",
                                      pressure_method = "Millero 1983")
    name == "automatic" && throw(ArgumentError(
        "K_method = \"automatic\" is no longer an option; use K_settings with " *
        RETIRED_K1_K2_NAMES["automatic"] * "."))
    return (K1_K2_method = name,)
end

"Whether a NamedTuple has any old-style method keyword. Settled from the type, so free."
@generated function _has_legacy(::NamedTuple{names}) where {names}
    return any(n in LEGACY_METHOD_KEYWORDS for n in names)
end

"""
    _fold_K_settings(arguments::NamedTuple) -> NamedTuple

Gather the constants' settings in `arguments` into a single `K_settings::KSettings`.

A NamedTuple `K_settings` becomes a `KSettings`. The old flat keywords (`K_method`,
`KF_method`, …) are translated, with a warning that shows the `K_settings` to write instead;
giving both forms is an error, since one would have to be ignored.

The translation of `K_method` wins over a legacy `BT_method`, as the old override did, so an
old call gets the numbers it always did.
"""
function _fold_K_settings(arguments::NamedTuple)
    if !_has_legacy(arguments)
        haskey(arguments, :K_settings) || return arguments
        arguments.K_settings isa KSettings && return arguments
        return merge(arguments, (K_settings = _as_K_settings(arguments.K_settings),))
    end

    legacy = NamedTuple{Tuple(k for k in keys(arguments) if k in LEGACY_METHOD_KEYWORDS)}(
        arguments)
    rest = Base.structdiff(arguments, legacy)
    haskey(rest, :K_settings) && throw(ArgumentError(
        "give the constants' settings either as K_settings or with the old keywords " *
        "($(join(keys(legacy), ", "))), not both"))

    flat = map(string, Base.structdiff(legacy, NamedTuple{(:K_method,)}))
    translated = haskey(legacy, :K_method) ? merge(flat, _translate_K_method(legacy.K_method)) :
                 flat
    settings = KSettings(; translated...)

    replacement = _non_default_fields(settings)
    @warn "$(join(keys(legacy), ", ")) $(length(legacy) == 1 ? "is" : "are") replaced by " *
          "K_settings. Write this instead:\n  K_settings = " *
          (isempty(replacement) ? "KSettings()" : _settings_literal(replacement)) *
          "\nThe old form still works for now." maxlog=1 _id=(:legacy_keywords, legacy)
    return merge(rest, (K_settings = settings,))
end

"""
What K1 and K2 must be multiplied by to move them from the scale they were fitted on to the
total scale. `sws_factor` converts seawater-scale to total; `fH` additionally relates NBS to
seawater-scale.
"""
_to_total_mult(scale::Symbol, sws_factor, fH) =
    scale === :sws ? sws_factor :
    scale === :nbs ? sws_factor / fH :
                     1.0    # :total — already there


"""
    calculate_constants(; temp_c, sal, pres_bar=0.0, K_settings=KSettings(), kwargs...)

Describe a body of seawater: its equilibrium constants, its totals, and the pH-scale
conversions between them.

Returns a NamedTuple `(; Ks, ST, FT, BT, Ca, Mg, fH, tf_fac, ts_fac, method)`, where `Ks` is
the inner bundle of 15 equilibrium constants and the rest is the water they were computed
for. Constants are on the total pH scale. `method` is the `KSettings` they were computed with.

```jldoctest
julia> env = CarbonateCalculator.calculate_constants(temp_c = 20.0, sal = 35.0);

julia> env.Ks.K1
1.2842425913487012e-6
```

Scalar only — the public entry points reject arrays before reaching here. To process many
samples, build a solver once and broadcast it; see [`CarbonateSystem`](@ref).

# Arguments

- `temp_c`: temperature, °C. Required.
- `sal`: practical salinity. Required.
- `pres_bar = 0.0`: hydrostatic pressure, bar.
- `ST`, `FT`, `BT = nothing`: total sulphate, fluoride and boron in mol/kg, derived from
  salinity when not given; `BT` follows `K_settings.BT_method`.
- `Ca`, `Mg = nothing`: calcium and magnesium in mol/kg. Modern seawater when not given; an
  explicit value overrides `K_settings.Ca_method` entirely.
- `K_settings = KSettings()`: which parameterisation to use for each constant, as a
  [`KSettings`](@ref) or a NamedTuple of its fields. See [`KSettings`](@ref) for every option.

The old keywords (`K_method`, `KSO4_method`, `BT_method`, `KF_method`, `KNH3_method`,
`Ca_method`, `MyAMI_mode`) are still accepted, with a warning naming the `K_settings` to use
instead.

# Extended help

The three parts — constants, totals, and scale-conversion factors — are returned together
because they are not independent. `tf_fac`, `ts_fac` and `fH` are built from the totals *and*
the constants, so a caller that recombined separately-derived pieces could convert a pH scale
using factors describing different water. Anyone who wants a single constant on its own
should call the parameterisation directly (`Lueker2000`, `calc_KB`, …) rather than take a
slice of this.

Every calculation in the package runs against one of these, and a solver given `Ks=` is given
*this* value rather than the inner bundle — see `_solve_core`.

Every choice comes from `K_settings` and from nothing else: no K1/K2 parameterisation
changes the totals, `fH` or the pressure correction it is used with. Under `"KGen"`, Kgen
supplies K0, K1, K2, KB, KW, KS, KF, the nutrient constants and both solubility products,
already pressure-corrected, so `KSO4_method` and `KF_method` have no effect and
`pressure_method` reaches only KH2S and KNH3.

See also [`KSettings`](@ref), [`CarbonateSystem`](@ref).
"""
function calculate_constants(; temp_c, sal, pres_bar=0.0, ST=nothing, FT=nothing,
    BT=nothing, Ca=nothing, Mg=nothing, K_settings=DEFAULT_K_SETTINGS, legacy...)

    settings = _settings_for_constants(K_settings, NamedTuple(legacy))
    K1_K2_method = settings.K1_K2_method

    # --- Boron, Sulfate, Fluoride and Calcium Totals ---
    # A supplied total always wins over the method that would otherwise derive it.
    final_BT = isnothing(BT) ? _apply(Val(:BT_method), settings.BT_method; sal).BT : BT
    final_ST = isnothing(ST) ? calc_ST(; sal).ST : ST
    final_FT = isnothing(FT) ? calc_FT(; sal).FT : FT

    # An unspecified Ca or Mg is modern seawater, and "modern" means the same composition
    # Kgen assumes for its KGen correction. Defaulting to that rather than to Culkin
    # keeps the calcium used for saturation state and the calcium used to correct the
    # constants the same number. Culkin and RT67 remain available explicitly.
    final_Ca = isnothing(Ca) ? _apply(Val(:Ca_method), settings.Ca_method; sal) : Ca
    final_Mg = something(Mg, MODERN_MAGNESIUM)

    # --- What every method needs -----------------------------------------------------------
    #
    # Kept out of the branches below because KGen's bundle, alone among the methods, arrives
    # from Kgen already complete and pressure-corrected. The only gaps in it are KH2S and
    # KNH3, which Kgen does not carry, so those two and their pressure factors are computed
    # for every method. Everything else the other parameterisations build on would be
    # computed and discarded on what is the default code path.
    pressure = _apply(Val(:pressure_method), settings.pressure_method)

    fH_val = _apply(Val(:fH_method), settings.fH_method; temp_c, sal).fH
    KH2S_surface = calc_KH2S(; temp_c, sal).KH2S
    KNH3_surface = _apply(Val(:KNH3_method), settings.KNH3_method; temp_c, sal).KNH3

    KH2S_at_pressure = KH2S_surface * _pressure_factor(temp_c, pres_bar, pressure.KH2S)
    KNH3_at_pressure = KNH3_surface * _pressure_factor(temp_c, pres_bar, pressure.KNH3)

    if K1_K2_method == "KGen"
        # Pass the totals computed above rather than deriving new ones here: those honour any
        # explicitly supplied ST/FT/BT and the choice of BT_method, and a second set of
        # formulas would discard BT_method on what is the default code path.

        # The seawater *composition* KGen corrects for: the explicit Ca/Mg if given, modern
        # otherwise. Deliberately not final_Ca, which may have come from Ca_method and be
        # salinity-scaled — that is a concentration, not a composition, and feeding it here
        # would read dilution as a low-Ca ocean.
        Mg_val = something(Mg, MODERN_MAGNESIUM)
        Ca_val = something(Ca, MODERN_CALCIUM)

        # Positional form, in Kgen's fixed order temp_c, sal, p_bar, magnesium, calcium.
        # Julia cannot broadcast keyword arguments, so this is the form that survives if this
        # branch is ever handed arrays directly. sulphate/fluorine/KGen_mode are keyword-only
        # in Kgen. Kgen owns the pressure correction, hence the pres_bar pass-through and the
        # absence of any correction block on this path.
        Ks = Kgen.calc_Ks(
            temp_c, sal, pres_bar, Mg_val, Ca_val;
            sulphate = final_ST,
            fluorine = final_FT,
            MyAMI_mode = settings.MyAMI_mode
        )

        # Kgen's bundle is on the total scale, so KNH3 needs the same conversion it gets on
        # the path below — built from Kgen's own KS and KF, which are the ones this bundle
        # describes. Skipping it leaves the default (Millero, seawater-scale) KNH3 about 2.2 %
        # high, visible only through NH4T, which defaults to zero.
        sws_factor = SWStoTOT(ST = final_ST, FT = final_FT, KS = Ks.KS, KF = Ks.KF)
        knh3_mult = _to_total_mult(_native_scale(Val(:KNH3_method), settings.KNH3_method),
                                   sws_factor, fH_val)

        final_Ks = merge(Ks, (KH2S = KH2S_at_pressure,
                              KNH3 = KNH3_at_pressure * knh3_mult))
    else
        # --- The surface environment the parameterisation is a correction to ---
        # `_surface` throughout: these are 1 atm values, not the ones that end up in the
        # bundle. K0 is the exception that never gets corrected, following CO2SYS.
        KS_surface = _apply(Val(:KSO4_method), settings.KSO4_method; temp_c, sal).KS
        KF_surface = _apply(Val(:KF_method), settings.KF_method; temp_c, sal).KF
        K0         = calc_K0(; temp_c, sal).K0
        KB_surface = calc_KB(; temp_c, sal, ST = final_ST, FT = final_FT,
                             KS = KS_surface, KF = KF_surface).KB
        KW_surface = calc_KW(; temp_c, sal).KW
        nutrients  = calc_KP(; temp_c, sal, fH = fH_val)

        # No parameterisation adjusts the environment above. Mehrbach 1973 and Millero 1979
        # used to replace parts of it (BT, fH, KB, KW, the nutrient constants and the
        # pressure table); those are now `KSettings` fields, chosen like any other, and
        # `KSettings` warns when the published companions are not selected.
        (; K1, K2) = _apply(Val(:K1_K2_method), K1_K2_method; temp_c, sal, ST = final_ST,
                            FT = final_FT, KS = KS_surface, KF = KF_surface, fH = fH_val)

        # --- Onto the total pH scale ---
        # KS and KF are pressure-corrected first, because the seawater-to-total factor is
        # built from the corrected values, not the surface ones.
        KS = KS_surface * _pressure_factor(temp_c, pres_bar, pressure.KS)
        KF = KF_surface * _pressure_factor(temp_c, pres_bar, pressure.KF)
        sws_factor = SWStoTOT(ST = final_ST, FT = final_FT, KS = KS, KF = KF)

        # KB and KH2S carry no multiplier because every parameterisation of them here is
        # total-scale native. K1/K2 and KNH3 read theirs off the table.
        k1k2_mult = _to_total_mult(_native_scale(Val(:K1_K2_method), K1_K2_method),
                                   sws_factor, fH_val)
        knh3_mult = _to_total_mult(_native_scale(Val(:KNH3_method), settings.KNH3_method),
                                   sws_factor, fH_val)

        final_Ks = (
            K1 = (K1 * _pressure_factor(temp_c, pres_bar, pressure.K1)) * k1k2_mult,
            K2 = (K2 * _pressure_factor(temp_c, pres_bar, pressure.K2)) * k1k2_mult,
            K0 = K0,
            KB = KB_surface * _pressure_factor(temp_c, pres_bar, pressure.KB),
            KS = KS,
            KF = KF,
            KH2S = KH2S_at_pressure,
            KW   = (KW_surface * _pressure_factor(temp_c, pres_bar, pressure.KW)) * sws_factor,
            KP1  = (nutrients.KP1 * _pressure_factor(temp_c, pres_bar, pressure.KP1)) * sws_factor,
            KP2  = (nutrients.KP2 * _pressure_factor(temp_c, pres_bar, pressure.KP2)) * sws_factor,
            KP3  = (nutrients.KP3 * _pressure_factor(temp_c, pres_bar, pressure.KP3)) * sws_factor,
            KSi  = (nutrients.KSi * _pressure_factor(temp_c, pres_bar, pressure.KSi)) * sws_factor,
            KNH3 = KNH3_at_pressure * knh3_mult,
            KspA = calc_KspA(temp_c = temp_c, sal = sal).KspA *
                   _pressure_factor(temp_c, pres_bar, pressure.KspA),
            KspC = calc_KspC(temp_c = temp_c, sal = sal).KspC *
                   _pressure_factor(temp_c, pres_bar, pressure.KspC)
        )
    end

    # pH scale conversion factors. These belong here beside `fH`, which is the same kind of
    # quantity for the NBS scale: they are pure functions of the totals and the constants
    # this function has just computed, and every caller needs them. Keeping them here is
    # what makes this return value a complete, self-consistent description of the water -
    # constants, totals and scale conversions together - so nothing downstream has to
    # recombine pieces of it and risk mixing supplied constants with derived totals.
    free_to_total = 1.0 + final_ST / final_Ks.KS
    sws_to_total = free_to_total / (1.0 + final_ST / final_Ks.KS + final_FT / final_Ks.KF)

    return (
        Ks = final_Ks,
        ST = final_ST,
        FT = final_FT,
        BT = final_BT,
        Ca = final_Ca,
        Mg = final_Mg,
        fH = fH_val,
        tf_fac = free_to_total,
        ts_fac = sws_to_total,
        method = settings
    )

end

"""
The `KSettings` a direct call to `calculate_constants` asked for.

The old flat keywords are folded the same way the solvers fold them, so a script that called
`calculate_constants(; K_method = "Lueker 2000")` keeps working, with the same warning. Anything
else in `legacy` is a typo, and an error: this is a public function and a misspelled setting
would otherwise be ignored.
"""
function _settings_for_constants(K_settings, legacy::NamedTuple)
    isempty(legacy) && return _as_K_settings(K_settings)
    unknown = [k for k in keys(legacy) if !(k in LEGACY_METHOD_KEYWORDS)]
    isempty(unknown) || throw(ArgumentError(
        "calculate_constants does not take $(join(unknown, ", "))"))
    given = K_settings === DEFAULT_K_SETTINGS ? legacy :
            merge(legacy, (K_settings = K_settings,))
    return _fold_K_settings(given).K_settings
end

export calculate_constants, KSettings
end # module