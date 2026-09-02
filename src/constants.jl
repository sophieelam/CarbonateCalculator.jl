"""
    Constants

Equilibrium constants, their temperature/salinity/pressure dependence, and the seawater
composition derived from salinity.

Every constant is available under several published parameterisations, selected by the
`*_method` arguments and registered in [`COMPOSITION`](@ref) and [`SPECIATION`](@ref).
[`calculate_constants`](@ref) assembles a complete set for one body of water;
[`which_K`](@ref) chooses a parameterisation when asked to.

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


# 16. Schockman & Byrne, 2020: https://doi.org/10.1016/j.gca.2021.02.008
# Salinity between 19.6-40 PSU, temperature between 15-35 C
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


# 17. Papdimitriou, et al., 2018: https://doi.org/10.1016/j.gca.2017.09.037
# Salinity between 50-100 PSU, temperature between (-6)-25 C
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
# Salinity between 19-43 PSU, temperature between 2-35 C
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
    BT=0
    return (; BT=BT) # in mol/kg-SW
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


# Case #8: Millero, 1979: https://doi.org/10.1016/0016-7037(79)90184-4
# Salinity of 0 PSU (freshwater), temperature between 0-50 C
function case8_fH(; temp_c, sal, kwargs...)
    fH = 1
    return (; fH=fH)
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

# Case #8: Millero, 1979: https://doi.org/10.1016/0016-7037(79)90184-4
# Salinity of 0 PSU (freshwater), temperature between 0-50 C
function case8_KB(; kwargs...)
    KB = 0
    return (; KB=KB)
end

# Case #6/7: Mehrbach 1973: https://doi.org/10.4319/lo.1973.18.6.0897
# Salinity between 19-43 PSU, temperature between 2-35 C
# All other cases
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

# Case #8: Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
# NOTES FROM CO2SYS:
# Refit data of Harned and Owen, The Physical Chemistry of Electrolyte Solutions, 1958
function case8_KW(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    lnKW = 148.9802 - 13847.26 / TK - 23.6521 * log(TK)
    KW = exp(lnKW)
    return (; KW=KW)
end



# AlL other cases:
# Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995
function calc_KW(; temp_c, sal, kwargs...)
    TK = temp_c + 273.15
    lnKW = (148.9802 - 13847.26 / TK - 23.6521 * log(TK) + (-5.977 + 118.67 /
    TK + 1.0495 * log(TK)) * sqrt(sal) - 0.01615 * sal)
    KW = exp(lnKW) # native SWS scale
    return (; KW=KW)
end


# Calculate KP1, KP2, KP3, and KSi

# Cases 6 & 8:
function case68_KP(; kwargs...)
    KP1 = 0
    KP2 = 0
    KP3 = 0
    KSi = 0
    # Neither the GEOSECS choice nor the freshwater choice
    # include contributions from phosphate or silicate.
    return (; KP1=KP1, KP2=KP2, KP3=KP3, KSi=KSi)
end


# All other cases
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


"""
    which_K(; K_method="default", temp_c, sal, Ca=nothing, Mg=nothing)

Pick the equilibrium-constant parameterisation best suited to the conditions.

Tested in two stages, and the order is important.

**Domain first.** `sal < 1` and `sal > 50` can only be modelled by Millero 1979 and
Papadimitriou 2018, which cover fresh water and brine respectively. Nothing may be preferred
ahead of them.

**Preference second.** Within the ordinary salinity range, the choice is between
parameterisations that all apply, so temperature decides.

Millero 2002 does take precedence over GP 1989 at `34 < sal < 37`. That is a genuine
preference between two parameterisations that both cover those conditions, not a domain error.
"""
function which_K(; K_method="default", temp_c, sal, Ca=nothing,
    Mg=nothing, kwargs...)
    # An unspecified Ca/Mg is modern seawater, so it must not read as "non-modern" here.
    if something(Ca, MODERN_CALCIUM) != MODERN_CALCIUM ||
       something(Mg, MODERN_MAGNESIUM) != MODERN_MAGNESIUM
        return "KGen"
    end

    # Domain guards.
    if sal < 1.0
        return "Millero 1979"
    elseif sal > 50.0
        return "Papadimitriou 2018"
    end

    # Preference guards, all within the salinity range these parameterisations share.
    if temp_c > 35
        return "Millero 2006"
    elseif temp_c < 2 && 34.0 < sal < 37.0
        return "Millero 2002"
    elseif temp_c < 0 && 10 < sal < 50
        return "GP 1989"
    else
        return "KGen"
    end
end

# --- Option lookup table ------------------------------------------------------------------

"""
Modern seawater calcium, as a parameterisation so that `Ca_method` has three interchangeable
options rather than two plus a fallback. Independent of salinity by definition: this is the
composition Kgen assumes for the KGen correction, not a concentration to be diluted.
"""
Modern_Ca(; kwargs...) = MODERN_CALCIUM

"""
How much of each dissolved species the water holds, by the name `*_method` selects.

Composition, not speciation: these parameterisations give a *concentration* from salinity.

Rows are `(fn = …)`. See [`SPECIATION`](@ref) for what the shared `default`/`options` shape
is for.
"""
const COMPOSITION = (
    BT_method = (
        default = "Uppstrom",
        options = Dict(
            "Uppstrom" => (fn = Uppstrom_BT,),
            "Lee" => (fn = Lee_BT,),
            "KSK18" => (fn = KSK18_BT,))
        ),
    Ca_method = (
        default = "modern",
        options = Dict(
            "modern" => (fn = Modern_Ca,),
            "Culkin" => (fn = Culkin_Ca,),
            "RT67" => (fn = RT67_Ca,))
        ),
)

"""
The equilibrium constants, by the name `*_method` selects, each with the pH scale it is
fitted on.

Rows are `(fn = …, scale = …)`

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
    # `default = nothing` because K_method's default is not a fixed name: `which_K` picks it
    # from temperature and salinity, and has already done so before `_apply` is reached.
    K_method    = (
        default = "KGen",
        options = Dict(
            "Roy 1993"           => (fn = Roy1993,           scale = :total),
            "GP 1989"            => (fn = GP1989,            scale = :sws),
            "Hansson 1973"       => (fn = Hansson1973,       scale = :sws),
            "DM 1987"            => (fn = DM1987,            scale = :sws),
            "HM 1973"            => (fn = HM1973,            scale = :sws),
            "Mehrbach 1973 A"    => (fn = Mehrbach1973,      scale = :sws),
            "Mehrbach 1973 B"    => (fn = Mehrbach1973,      scale = :sws),
            "Millero 1979"       => (fn = Millero1979,       scale = :nbs),
            "CW 2003"            => (fn = CW2003,            scale = :sws),
            "Lueker 2000"        => (fn = Lueker2000,        scale = :total),
            "MPM 2002"           => (fn = MPM2002,           scale = :sws),
            "Millero 2002"       => (fn = Millero2002,       scale = :sws),
            "Millero 2006"       => (fn = Millero2006,       scale = :sws),
            "Millero 2010"       => (fn = Millero2010,       scale = :sws),
            "Waters 2014"        => (fn = Waters2014,        scale = :sws),
            "SB 2020"            => (fn = SB2020,            scale = :total),
            "Papadimitriou 2018" => (fn = Papadimitriou2018, scale = :total),
            "Sulpis 2020"        => (fn = Sulpis2020,        scale = :total),
            # Listed so this is the complete set of valid K_method values and the error message
            # stays right, but it has no K1/K2 function: Kgen hands back a whole bundle, already
            # pressure-corrected and on the total scale, assembled separately below.
            "KGen"              => (fn = nothing,           scale = :total),
        )
    ),
)

"""
Raise for a `*_method` value that names no parameterisation.

`valid` is built from the table when the dispatch below is generated, so the list of names
cannot drift from the ones actually accepted, and the error path touches no table at all.
"""
_unknown_method(argument::Symbol, name, valid::String) =
    throw(ArgumentError("Unknown $argument: \"$name\".\nChoose one of: $valid"))

"The 'you could have meant' list for one argument, built at generation time."
function _valid_names(entry)
    fallback = isnothing(entry.default) ?
        "\"default\" picks one from temperature and salinity" :
        "\"default\" is $(entry.default)"
    return join(sort!(collect(keys(entry.options))), ", ") * ". $fallback."
end

"""
    _apply(Val(:KSO4_method), name; kwargs...)

Call the parameterisation `name` selects for one `*_method` argument, resolving `"default"`
and rejecting anything the table does not list.

**Generated from [`COMPOSITION`](@ref) and [`SPECIATION`](@ref), not a lookup into them.** A
`Dict` holding several different functions can only be typed `Dict{String,Function}`, so
calling through it is a dynamic dispatch returning `Any` — and that `Any` spreads into every
constant derived from the result, leaving the whole bundle untyped. Interpolating the function
objects into a comparison chain at load time keeps every call statically known while the
tables stay the single place a parameterisation is registered. Do not replace this with a
`Dict` lookup - it is much slower.
"""
function _apply end

"""
    _native_scale(Val(:K_method), name) -> :total | :sws | :nbs

The pH scale the parameterisation `name` selects is fitted on.

Generated from [`SPECIATION`](@ref) only — the composition parameterisations give
concentrations, which are not on a pH scale.

Resolves `"default"` exactly as [`_apply`](@ref) does, which is the point of going through it
at all: a scale test written against the raw argument, such as
`KNH3_method == "Millero" ? sws_factor : 1.0`, silently misses the `"default"` that selects
Millero, and leaves a seawater-scale constant in a total-scale bundle.
"""
function _native_scale end

for table in (COMPOSITION, SPECIATION), (argument, entry) in pairs(table)
    body = Expr(:block)

    # A default that is a fixed name resolves here. K_method's is not — `which_K` derives it
    # from temperature and salinity long before this is reached — so it emits no such line.
    isnothing(entry.default) ||
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
    isnothing(entry.default) ||
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
What K1 and K2 must be multiplied by to move them from the scale they were fitted on to the
total scale. `sws_factor` converts seawater-scale to total; `fH` additionally relates NBS to
seawater-scale.
"""
_to_total_mult(scale::Symbol, sws_factor, fH) =
    scale === :sws ? sws_factor :
    scale === :nbs ? sws_factor / fH :
                     1.0    # :total — already there


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
Millero 1983's coefficients, which the freshwater `Millero 1979` parameterisation uses in
place of the rows above.

A table rather than a branch, so the special case is three rows that `merge` overlays on
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

"""
The pressure coefficients to use for `K_method`.

Identical field names and types either way, so the choice stays inferable — a `Union` of two
differently-shaped tables here would leave every corrected constant untyped.
"""
_pressure_table(K_method) = K_method == "Millero 1979" ?
    merge(PRESSURE_COEFFICIENTS, MILLERO1979_PRESSURE) : PRESSURE_COEFFICIENTS

"""
    calculate_constants(; temp_c, sal, pres_bar=0.0, kwargs...)

Describe a body of seawater: its equilibrium constants, its totals, and the pH-scale
conversions between them.

Returns a NamedTuple `(; Ks, ST, FT, BT, Ca, Mg, fH, tf_fac, ts_fac, method)`, where `Ks` is
the inner bundle of 15 equilibrium constants and the rest is the water they were computed
for. Constants are on the total pH scale.

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
  salinity when not given.
- `Ca`, `Mg = nothing`: calcium and magnesium in mol/kg. Modern seawater when not given; an
  explicit value overrides `Ca_method` entirely.
- `MyAMI_mode = "approximate"`: how the KGen path corrects the constants for Ca and Mg.
  `"calculate"` runs the full MyAMI model, which is slower and more accurate.

Each `*_method` selects a parameterisation, and `"default"` picks the recommended one.

`K_method` selects the carbonic-acid parameterisation, and defaults to `"KGen"`. Its options are
`"KGen"`, `"Roy 1993"`, `"GP 1989"`, `"Hansson 1973"`, `"DM 1987"`, `"HM 1973"`,
`"Mehrbach 1973 A"`, `"Mehrbach 1973 B"`, `"Millero 1979"`, `"CW 2003"`, `"Lueker 2000"`,
`"MPM 2002"`, `"Millero 2002"`, `"Millero 2006"`, `"Millero 2010"`, `"Waters 2014"`,
`"SB 2020"`, `"Papadimitriou 2018"`, `"Sulpis 2020"` and `"automatic"`, which picks the one
best suited to the conditions — see [`which_K`](@ref).

The rest:

| argument | options | default |
|---|---|---|
| `KSO4_method` | `"Dickson"`, `"Khoo"`, `"WM13"` | `"Dickson"` (CO2SYS's recommendation) |
| `BT_method` | `"Uppstrom"`, `"Lee"`, `"KSK18"` | `"Uppstrom"` (CO2SYS's recommendation) |
| `KF_method` | `"Dickson"`, `"Perez"` | `"Dickson"` (CO2SYS's recommendation) |
| `KNH3_method` | `"Millero"`, `"Clegg"` | `"Millero"` |
| `Ca_method` | `"Culkin"`, `"RT67"` | modern seawater Ca (`MODERN_CALCIUM`) |

`Ca_method`'s default is modern seawater rather than Culkin because that is the composition
Kgen assumes for the MyAMI correction, so the Ca used for saturation state and the Ca used to
correct the constants agree. `"Culkin"` and `"RT67"` give salinity-scaled concentrations
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

Pressure corrections follow Millero (1995), except under `"Millero 1979"`, which uses the
1983 freshwater coefficients, and under `"KGen"`, where Kgen owns the correction.

See also [`which_K`](@ref), [`CarbonateSystem`](@ref).
"""
function calculate_constants(; temp_c, sal, pres_bar=0.0, ST=nothing, FT=nothing,
    BT=nothing, K_method="KGen", KSO4_method="Dickson", BT_method="Uppstrom",
    KF_method="Dickson", KNH3_method="Millero", Ca_method="modern", kwargs...)

    # `"default"` is resolved here rather than in `_apply`, because the method it names has no
    # K1/K2 function to dispatch to — Kgen delivers a whole bundle, and the branch below is what
    # reaches it. Left to `_apply`, "default" would be rewritten to "KGen" and then fail to
    # match any clause, which reports an unknown method that is in its own list of valid names.
    if K_method == "default"
        K_method = SPECIATION.K_method.default
    elseif K_method == "automatic"
        K_method = which_K(K_method=K_method, temp_c=temp_c, sal=sal; kwargs...)
    end

    # --- Boron, Sulfate, Fluoride and Calcium Totals ---
    # A supplied total always wins over the method that would otherwise derive it.
    final_BT = isnothing(BT) ? _apply(Val(:BT_method), BT_method; sal).BT : BT
    final_ST = isnothing(ST) ? calc_ST(; sal).ST : ST
    final_FT = isnothing(FT) ? calc_FT(; sal).FT : FT

    Ca_in = get(kwargs, :Ca, nothing)
    Mg_in = get(kwargs, :Mg, nothing)

    # An unspecified Ca or Mg is modern seawater, and "modern" means the same composition
    # Kgen assumes for its KGen correction. Defaulting to that rather than to Culkin
    # keeps the calcium used for saturation state and the calcium used to correct the
    # constants the same number. Culkin and RT67 remain available explicitly.
    final_Ca = isnothing(Ca_in) ? _apply(Val(:Ca_method), Ca_method; sal) : Ca_in
    final_Mg = something(Mg_in, MODERN_MAGNESIUM)

    # --- What every method needs -----------------------------------------------------------
    #
    # Kept out of the branches below because KGen's bundle, alone among the methods, arrives
    # from Kgen already complete and pressure-corrected. The only gaps in it are KH2S and
    # KNH3, which Kgen does not carry, so those two and their pressure factors are computed
    # for every method. Everything else the other parameterisations build on would be
    # computed and discarded on what is the default code path.
    pressure = _pressure_table(K_method)

    fH_val = calc_fH(temp_c, sal)
    KH2S_surface = calc_KH2S(; temp_c, sal).KH2S
    KNH3_surface = _apply(Val(:KNH3_method), KNH3_method; temp_c, sal).KNH3

    KH2S_at_pressure = KH2S_surface * _pressure_factor(temp_c, pres_bar, pressure.KH2S)
    KNH3_at_pressure = KNH3_surface * _pressure_factor(temp_c, pres_bar, pressure.KNH3)

    if K_method == "KGen"
        # Pass the totals computed above rather than deriving new ones here: those honour any
        # explicitly supplied ST/FT/BT and the choice of BT_method, and a second set of
        # formulas would discard BT_method on what is the default code path.

        # The seawater *composition* KGen corrects for: the explicit Ca/Mg if given, modern
        # otherwise. Deliberately not final_Ca, which may have come from Ca_method and be
        # salinity-scaled — that is a concentration, not a composition, and feeding it here
        # would read dilution as a low-Ca ocean.
        Mg_val = something(Mg_in, MODERN_MAGNESIUM)
        Ca_val = something(Ca_in, MODERN_CALCIUM)

        MyAMI_mode = get(kwargs, :MyAMI_mode, "approximate")

        # Positional form, in Kgen's fixed order temp_c, sal, p_bar, magnesium, calcium.
        # Julia cannot broadcast keyword arguments, so this is the form that survives if this
        # branch is ever handed arrays directly. sulphate/fluorine/KGen_mode are keyword-only
        # in Kgen. Kgen owns the pressure correction, hence the pres_bar pass-through and the
        # absence of any correction block on this path.
        Ks = Kgen.calc_Ks(
            temp_c, sal, pres_bar, Mg_val, Ca_val;
            sulphate = final_ST,
            fluorine = final_FT,
            MyAMI_mode = MyAMI_mode
        )

        # Kgen's bundle is on the total scale, so KNH3 needs the same conversion it gets on
        # the path below — built from Kgen's own KS and KF, which are the ones this bundle
        # describes. Skipping it leaves the default (Millero, seawater-scale) KNH3 about 2.2 %
        # high, visible only through NH4T, which defaults to zero.
        sws_factor = SWStoTOT(ST = final_ST, FT = final_FT, KS = Ks.KS, KF = Ks.KF)
        knh3_mult = _to_total_mult(_native_scale(Val(:KNH3_method), KNH3_method),
                                   sws_factor, fH_val)

        final_Ks = merge(Ks, (KH2S = KH2S_at_pressure,
                              KNH3 = KNH3_at_pressure * knh3_mult))
    else
        # --- The surface environment the parameterisation is a correction to ---
        # `_surface` throughout: these are 1 atm values, not the ones that end up in the
        # bundle. K0 is the exception that never gets corrected, following CO2SYS.
        KS_surface = _apply(Val(:KSO4_method), KSO4_method; temp_c, sal).KS
        KF_surface = _apply(Val(:KF_method), KF_method; temp_c, sal).KF
        K0         = calc_K0(; temp_c, sal).K0
        KB_surface = calc_KB(; temp_c, sal, ST = final_ST, FT = final_FT,
                             KS = KS_surface, KF = KF_surface).KB
        KW_surface = calc_KW(; temp_c, sal).KW
        nutrients  = calc_KP(; temp_c, sal, fH = fH_val)

        # Three parameterisations refit part of that surface environment before their K1/K2
        # can be used. They are spelled out because each replaces a different set; the other
        # fifteen need nothing but the lookup. Note that `nutrients` above keeps the baseline
        # fH even where fH is refitted here — only Millero 1979 replaces those constants.
        if K_method == "Mehrbach 1973 A"
            final_BT = isnothing(BT) ? case67_BT(; sal).BT : BT
        elseif K_method == "Mehrbach 1973 B"
            fH_val = case7_fH(; temp_c, sal).fH
        elseif K_method == "Millero 1979"
            # Freshwater, so boron and the nutrient constants go to zero rather than being
            # rescaled.
            fH_val = case8_fH(; temp_c, sal).fH
            final_BT = isnothing(BT) ? case8_BT(; sal).BT : BT
            KB_surface = case8_KB(; temp_c, sal, fH_val).KB
            KW_surface = case8_KW(; temp_c, sal).KW
            nutrients = case68_KP(; temp_c, sal, fH_val)
        end

        (; K1, K2) = _apply(Val(:K_method), K_method; temp_c, sal, ST = final_ST,
                            FT = final_FT, KS = KS_surface, KF = KF_surface, fH = fH_val)

        # --- Onto the total pH scale ---
        # KS and KF are pressure-corrected first, because the seawater-to-total factor is
        # built from the corrected values, not the surface ones.
        KS = KS_surface * _pressure_factor(temp_c, pres_bar, pressure.KS)
        KF = KF_surface * _pressure_factor(temp_c, pres_bar, pressure.KF)
        sws_factor = SWStoTOT(ST = final_ST, FT = final_FT, KS = KS, KF = KF)

        # KB and KH2S carry no multiplier because every parameterisation of them here is
        # total-scale native. K1/K2 and KNH3 read theirs off the table.
        k1k2_mult = _to_total_mult(_native_scale(Val(:K_method), K_method), sws_factor, fH_val)
        knh3_mult = _to_total_mult(_native_scale(Val(:KNH3_method), KNH3_method),
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
        method = K_method
    )

end

export calculate_constants
end # module