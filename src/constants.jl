module Constants
using PythonCall

# 1. Create empty placeholders instead of importing right away
const np = PythonCall.pynew()
const kgen = PythonCall.pynew()

# 2. Fill them exactly when the module is loaded at runtime
function __init__()
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(kgen, pyimport("kgen"))
end

using ..Helpers
using Statistics


# Helper function:
function SWStoTOT(; ST, FT, KS, KF, kwargs...)
    return (1 + ST / KS) / (1 + ST / KS + FT / KF)
end

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
function Roy1993(; T_in, S_in, ST, FT, KS, KF, kwargs...)
    TK = T_in + 273.15
    logTK = log(TK)
    sqrS = sqrt(S_in)
    S15 = S_in * sqrS

    # K1 Calculation
    lnK1 = (2.83655 - 2307.1266 / TK - 1.5529413 * logTK +
           (-0.20760841 - 4.0484 / TK) * sqrS + 0.08468345 * S_in -
           0.00654208 * S15)
    
    # Originally on total scale! "* (1 - 0.001005109 * S_in)" converts to mol/kg-sw
    # then, if I were to "/ SWStoTOT" I would get to SWS scale
    K1 = exp(lnK1) * (1 - 0.001005109 * S_in)

    # K2 Calculation 
    lnK2 = (-9.226508 - 3351.6106 / TK - 0.2005743 * logTK +
           (-0.106901773 - 23.9722 / TK) * sqrS + 0.1130822 * S_in -
           0.00846934 * S15)
    
    # Originally on total scale! "* (1 - 0.001005109 * S_in)" converts to mol/kg-sw
    # then, if I were to "/ SWStoTOT" I would get to SWS scale
    K2 = exp(lnK2) * (1 - 0.001005109 * S_in)

    return (; K1=K1, K2=K2)
end


# 2. Goyet & Poisson, 1989: https://doi.org/10.1016/0198-0149(89)90064-2
# Salinity between 10-50 PSU, temperature between -1-40 C
# Seawater pH scale, artificial seawater
# NOTES FROM CO2SYS:
# The 2s precision in pK1 is .011, or 2.5% in K1.
# The 2s precision in pK2 is .02, or 4.5% in K2.
# This is in Table 5 on p. 1652 and what they use in the abstract:
function GP1989(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = (812.27 / TK + 3.356 - 0.00171 * S_in * log(TK) +
    0.000091 * S_in^2)

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = (1450.87 /TK + 4.604 - 0.00385 * S_in * log(TK) +
    0.000182 * S_in^2)

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
function Hansson1973(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = 851.4 /TK + 3.237 - 0.0106 * S_in + 0.000105 * S_in^2

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = (-3885.4 / TK + 125.844 - 18.141 * log(TK) - 0.0192 * S_in +
    0.000132 * S_in^2)

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
function DM1987(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = (3670.7 / TK - 62.008 + 9.7944 * log(TK) - 0.0118 * S_in +
    0.000116 * S_in^2)

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = 1394.7 / TK + 4.777 - 0.0184 * S_in + 0.000118 * S_in^2

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
function HM1973(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = 845 / TK + 3.248 - 0.0098 * S_in + 0.000087 * S_in^2

    # On the SWS pH scale in mol/kg-sw
    K1 = 10.0^(-pK1)

    pK2 = 1377.3 / TK + 4.824 - 0.0185 * S_in + 0.000122 * S_in^2

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
function Mehrbach1973(; T_in, S_in, fH, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = (-13.7201 + 0.031334 * TK + 3235.76 / TK + 1.3e-5 * S_in *
    TK - 0.1032 * sqrt(S_in))

    # On the SWS pH scale in mol/kg-sw (converted from NBS scale)
    K1 = 10.0^(-pK1) / fH

    pK2 = (5371.9645 + 1.671221 * TK + 0.22913 * S_in + 18.3802 *
    log10(S_in)) - 128375.28 / TK - 2194.3055 * log10(TK) -
    8.0944e-4 * S_in * TK - 5617.11 * log10(S_in) / TK + 2.136 *
    S_in / TK

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
function Millero1979(; T_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

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
function CW2003(; T_in, S_in, fH, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin
    F1 = 200.1 / TK + 0.3220
    F2 = -129.24 / TK + 1.4381

    pK1 = 3404.71 / TK + 0.032786 * TK - 14.8435 - 0.071692 * F1 * 
    sqrt(S_in) + 0.0021487 * S_in

    K1 = 10.0^(-pK1) / fH

    pK2 =  2902.39 / TK + 0.02379 * TK - 6.4980 - 0.3191 * F2 * 
    sqrt(S_in) + 0.0198 * S_in

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
function Lueker2000(; T_in, S_in, ST, FT, KS, KF, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = 3633.86 / TK - 61.2172 + 9.6777 * log(TK) - 0.011555 *
    S_in + 0.0001152 * S_in^2

    K1 = 10.0^(-pK1) # / SWStoTOT(ST=ST, FT=FT, KS=KS, KF=KF) # converted to SWS

    pK2 =  471.78 / TK + 25.929 - 3.16967 * log(TK) - 0.01781 *
    S_in + 0.0001122 * S_in^2

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
function MPM2002(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = -43.6977 - 0.0129037 * S_in + 1.364e-4 * S_in^2 +
    2885.378 / TK + 7.045159 * log(TK)

    K1 = 10.0^(-pK1) 

    pK2 = (-452.0940 + 13.142162 * S_in - 8.101e-4 * S_in^2 +
    21263.61 / TK + 68.483143 * log(TK) + (-581.4428 * S_in +
    0.259601 * S_in^2) / TK - 1.967035 * S_in * log(TK))

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
function Millero2002(; T_in, S_in, kwargs...)

    pK1 = 6.359 - 0.00664 * S_in - 0.01322 * T_in + 4.989e-5 * T_in^2

    K1 = 10.0^(-pK1) 

    pK2 = 9.867 - 0.01314 * S_in - 0.01904 * T_in + 2.448e-5 * T_in^2

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end 


# 13. Millero, et al., 2006: https://doi.org/10.1016/j.marchem.2005.12.001
# Salinity between 1-50 PSU, temperature between 0-50 C
# pH on seawater scale, real seawater
function Millero2006(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    pK2₀ = -90.18333 + 5143.692 / TK + 14.613358 * log(TK)
    A₁ = 13.4191 * sqrt(S_in) + 0.0331 * S_in - 5.33e-5 * S_in^2
    B₁ = -530.123 * sqrt(S_in) - 6.103 * S_in
    C₁ = -2.06950 * sqrt(S_in)
    A₂ = 21.0894 * sqrt(S_in) + 0.1248 * S_in - 3.687e-4 * S_in^2
    B₂ = -772.483 * sqrt(S_in) - 20.051 * S_in
    C₂ = -3.3336 * sqrt(S_in)

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
function Millero2010(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    pK2₀ = -90.18333 + 5143.692 / TK + 14.613358 * log(TK)
    A₁ = 13.4038 * sqrt(S_in) + 0.03206 * S_in - 5.242e-5 * S_in^2
    B₁ = -530.659 * sqrt(S_in) - 5.8210 * S_in
    C₁ = -2.0664 * sqrt(S_in)
    A₂ = 21.3728 * sqrt(S_in) + 0.1218 * S_in - 3.688e-4 * S_in^2
    B₂ = -788.289 * sqrt(S_in) - 19.189 * S_in
    C₂ = -3.374 * sqrt(S_in)

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
function Waters2014(; T_in, S_in, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    pK2₀ = -90.18333 + 5143.692 / TK + 14.613358 * log(TK)
    A₁ = 13.409160 * sqrt(S_in) + 0.031646 * S_in - 5.1895e-5 * S_in^2
    B₁ = -531.3642 * sqrt(S_in) - 5.713 * S_in
    C₁ = -2.0669166 * sqrt(S_in)
    A₂ = 21.225890 * sqrt(S_in) + 0.12450870 * S_in - 3.7243e-4 * S_in^2
    B₂ = -779.3444 * sqrt(S_in) - 19.91739 * S_in
    C₂ = -3.3534679 * sqrt(S_in)

    pK1 = pK1₀ + A₁ + B₁ / TK + C₁ * log(TK)

    K1 = 10.0^(-pK1)

    pK2 = pK2₀ + A₂ + B₂ / TK + C₂ * log(TK)

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end


# 16. Schockman & Byrne, 2020: https://doi.org/10.1016/j.gca.2021.02.008
# Salinity between 19.6-40 PSU, temperature between 15-35 C
# pH on total scale, converted to SWS if /SWStoTOT is included
function SB2020(; T_in, S_in, ST=nothing, FT=nothing, KS=nothing, KF=nothing, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin
    
    # K1 Reference: Waters, Millero, and Woosley (2014) coefficients used by MATLAB
    pK1₀ = -126.34048 + 6320.813 / TK + 19.568224 * log(TK)
    A₁   = 13.568513 * sqrt(S_in) + 0.031645 * S_in - 5.3834e-5 * S_in^2
    B₁   = -539.2304 * sqrt(S_in) - 5.635 * S_in
    C₁   = -2.0901396 * sqrt(S_in)

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

    pK2 = e1 + e2 / TK + e3 * log(TK) + e4 * S_in + e5 * sqrt(S_in) +
          e6 * S_in^2 + e7 * (S_in / TK)

    K2 = 10.0^(-pK2) 

    return (; K1=K1, K2=K2)
end


# 17. Papdimitriou, et al., 2018: https://doi.org/10.1016/j.gca.2017.09.037
# Salinity between 50-100 PSU, temperature between (-6)-25 C
# pH on total scale, converted to SWS if converted to SWS if /SWStoTOT is included
function Papadimitriou2018(; T_in, S_in, ST, FT, KS, KF, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = -176.48 + 6.14528 * sqrt.(S_in) - 0.127714 * S_in + 7.396 * 
    10.0^(-5) * S_in^2 + (9914.37 - 622.886 * sqrt(S_in) + 29.714 * S_in) /
    TK + (26.05129 - 0.666812 * sqrt(S_in)) * log(TK)

    pK2 = -323.52692 + 27.557655 * sqrt(S_in) + 0.154922 * S_in -
    2.48396 * 10.0^(-4) * S_in^2 + (14763.287 - 1014.819 * sqrt(S_in) -
    14.35223 * S_in) / TK + (50.385807 - 4.4630415 * sqrt(S_in)) * log(TK)

    K1 = 10.0^(-pK1) # / SWStoTOT(ST = ST, FT=FT, KS=KS, KF=KF) 

    K2 = 10.0^(-pK2) # / SWStoTOT(ST = ST, FT=FT, KS=KS, KF=KF) # converted to SWS 

    return (; K1=K1, K2=K2)
end


# 18. Sulpis, et al., 2020: https://doi.org/10.5194/os-16-847-2020
# Salinity between 19-43 PSU, temperature between 2-35 C
# Total pH scale, real seawater (converted to SWS converted to SWS if /SWStoTOT is included)
function Sulpis2020(; T_in, S_in, ST, FT, KS, KF, kwargs...)
    TK = T_in + 273.15     # Convert to Kelvin

    pK1 = 8510.63 / TK - 172.4493 + 26.32996 * log(TK) - 0.011555 * 
    S_in + 0.0001152 * S_in^2

    K1 = 10.0^(-pK1) # / SWStoTOT(ST=ST, FT=FT, KS=KS, KF=KF)

    pK2 =  4226.23 / TK - 59.4636 + 9.60817 * log(TK) - 0.01781 *
    S_in + 0.0001122 * S_in^2

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
function case67_BT(; S_in, kwargs...)
    BT = 0.0004106 * S_in / 35.0 # in mol/kg-SW

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
function Uppstrom_BT(; S_in, kwargs...)
    BT = 0.0004157 * S_in / 35.0 # in mol/kg-SW
    return (; BT=BT) # in mol/kg-SW
end

# Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.
# Geochimica Et Cosmochimica Acta 74 (6): 1801â€“1811.
function Lee_BT(; S_in, kwargs...)
    BT = 0.0004326 * S_in / 35.0 # in mol/kg-SW

    return (; BT=BT) # in mol/kg-SW
end

# NOTES FROM CO2SYS:
# Note that the reference provided an equation for μmol/kg-sw
# This function divides it by a factor of 1e6 to convert to mol/kg-sw
function KSK18_BT(; S_in, kwargs...)
    return (10.838 * S_in + 13.821) / 1e6 # in mol/kg-SW
end

# Calculate FT
# Riley, J. P., Deep-Sea Research 12:219-220, 1965
function  calc_FT(; S_in, kwargs...)
    FT = (0.000067 / 18.998) * (S_in / 1.80655)
    return (; FT=FT) # in mol/kg-SW
end

# Calculate ST
# Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966
function calc_ST(; S_in, kwargs...)
    ST = (0.14 / 96.062) * (S_in / 1.80655)
    return (; ST=ST) # in mol/kg-SW
end


# Calculate K0
# Weiss, R. F., Marine Chemistry 2:203-215, 1974
function calc_K0(; T_in, S_in, kwargs...)
    TempK100 = (T_in + 273.15) / 100
    lnK0 = (-60.2409 + 93.4517 / TempK100 + 23.3585 * log(TempK100) + 
    S_in * (0.023517 - 0.023656 * TempK100 + 0.0047036 * TempK100^2))

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

function Dickson_KS(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    # Ionic strength on the molal scale
    IonS = 19.924 * S_in / (1000.0 - 1.005 * S_in)

    lnKS = (-4276.1 / T + 141.328 - 23.093 * log(T)) + 
           (-13856.0 / T + 324.57 - 47.986 * log(T)) * sqrt(IonS) + 
           (35474.0 / T - 771.54 + 114.723 * log(T)) * IonS + 
           (-2698.0 / T) * IonS^1.5 + (1776.0 / T) * IonS^2

    # Convert from mol/kg-H2O to mol/kg-SW
    KS = exp(lnKS) * (1.0 - 0.001005109 * S_in) 
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
function Khoo_KS(; T_in, S_in, kwargs...)
    T = T_in + 273.15

    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
    IonS = 19.924 * S_in / (1000 - 1.005 * S_in)

    pKS = 647.59 / T - 6.3451 + 0.019085 * T - 0.5208 * sqrt(IonS)
    KS = (10.0^(-pKS)) * (1 - 0.001005109 * S_in) # mol/kg-SW

    return (; KS=KS)

end

function WM13_KS(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    logKS0 = (562.69486 - 102.5154 * log(T) - 0.0001117033 * T^2 + 0.2477538 * T -
    13273.76 / T)

    logKSK0 = ((4.24666 - 0.152671 * T + 0.0267059 * T * log(T) - 0.000042128 *
    T^2) * S_in^0.5 + (0.2542181 - 0.00509534 * T + 0.00071589 * T * log(T) * S_in) +
    (-0.00291179 + 0.0000209968 * T) * S_in^1.5 -0.0000403724 * S_in^2)

    kSO4 = (1 - 0.001005109 * S_in) * 10^(logKS0 + logKSK0)
    return kSO4 # in mol/kg-SW
end

# Calculate KF

# Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979
function Dickson_KF(; T_in, S_in, kwargs...)
    T = T_in + 273.15

    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
    IonS = 19.924 * S_in / (1000 - 1.005 * S_in)
    lnKF = 1590.2 / T - 12.641 + 1.525 * sqrt(IonS)
    KF = (exp(lnKF) * (1.0 - 0.001005109 * S_in)) # in mol/kg-SW
    return (; KF=KF)
end

# Perez & Fraga, 1987: https://doi.org/10.1016/0304-4203(87)90036-3
# Salinity 10-40 PSU and Temperature of 9-33 C
# NOTE FrOM PYCO2SYS:
# Note that this is not currently used or an option in CO2SYS,
# despite the equations below appearing in CO2SYS.m (commented out).
function Perez_KF(; T_in, S_in, kwargs...)
    T = T_in + 273.15

    lnKF = 874 / T - 9.68 + 0.111 * sqrt(S_in)
    KF = (exp(lnKF) * (1.0 - 0.001005109 * S_in)) # # in mol/kg-SW
    return (; KF=KF)
end



# Calculate fH

# Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.

# Case #8: Millero, 1979: https://doi.org/10.1016/0016-7037(79)90184-4
# Salinity of 0 PSU (freshwater), temperature between 0-50 C
function case8_fH(; T_in, S_in, kwargs...)
    fH = 1
    return (; fH=fH)
end

# Case #7: Mehrbach 1973: https://doi.org/10.4319/lo.1973.18.6.0897
# Salinity between 19-43 PSU, temperature between 2-35 C
# NBS pH scale converted to SWS, real seawater
function case7_fH(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    fH = (1.29 - 0.00204 * T + (0.00046 - 0.00000148 * T) * S_in * S_in)
    
    # NOTES FROM CO2SYS
    # Peng et al, Tellus 39B:439-458, 1987:
    # They reference the GEOSECS report, but round the value
    # given there off so that it is about .008 (1#) lower. It
    # doesn't agree with the check value they give on p. 456.
    return (; fH=fH)
end


# All other cases
# Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80)
function calc_fH(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    fH = (1.2948 - 0.002036 * T + (0.0004607 - 0.000001475 * T) * S_in^2)
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
# NBS pH scale converted to SWS, real seawater
# NOTES FROM CO2SYS:
# This is for GEOSECS and Peng et al.
# Lyman, John, UCLA Thesis, 1957 fit by Li et al, JGR 74:5507-5525, 1969:
function case67_KB(; T_in, S_in, fH, kwargs...)
    while fH isa NamedTuple
        fH = fH
    end
    logKB = -9.26 + 0.00886 * S_in + 0.01 * T_in
    KB = (10.0^(logKB) / fH) # converted to SWS pH scale
    return (; KB=KB)
end


# All other cases
# Dickson, A. G., Deep-Sea Research 37:755-766, 1990
function calc_KB(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    lnKBnum = (-8966.9 - 2890.53 * sqrt(S_in) - 77.942 * S_in + 1.728 *
    sqrt(S_in) * S_in - 0.0996 * S_in^2)
    lnKB = (lnKBnum / T + 148.0248 + 137.1942 * sqrt(S_in) + 1.62142 *
    S_in + (-24.4344 - 25.085 * sqrt(S_in) - 0.2474 * S_in) * log(T) +
    0.053105 * sqrt(S_in) * T)

    KB = (exp(lnKB)) # native total scale
    return (; KB=KB)
end


# Calculate KW

# Case #7: Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
function case7_KW(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    lnKW = (148.9802 - 13847.26 / T - 23.6521 * log(T) + (-79.2447 + 
    3298.72 / T + 12.0408 * log(T)) * sqrt(S_in) - 0.019813 * S_in)
    KW = exp(lnKW)
    return (; KW=KW)
end


# Case #8: Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
# NOTES FROM CO2SYS:
# Refit data of Harned and Owen, The Physical Chemistry of Electrolyte Solutions, 1958
function case8_KW(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    lnKW = 148.9802 - 13847.26 / T - 23.6521 * log(T)
    KW = exp(lnKW)
    return (; KW=KW)
end



# Case #6: 
function case6_KW(; T_in, S_in, kwargs...)
    KW = 0 # GEOSECS doesn't include OH effects
    return (; KW=KW)
end


# AlL other cases:
# Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995
function calc_KW(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    lnKW = (148.9802 - 13847.26 / T - 23.6521 * log(T) + (-5.977 + 118.67 /
    T + 1.0495 * log(T)) * sqrt(S_in) - 0.01615 * S_in)
    KW = exp(lnKW) # native SWS scale
    return (; KW=KW)
end


# Calculate KP1, KP2, KP3, and KSi

# Case #7: Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
function case7_KP(; T_in, S_in, fH, kwargs...)
    T = T_in + 273.15
    KP1 = 0.02
    # NOTE FROM CO2SYS:
    # Peng et al don't include the contribution from this term,
    # but it is so small it doesn't contribute. It needs to be
    # kept so that the routines work ok.
    # KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
    # Limnology and Oceanography 12:243-252, 1967
    # these are only for salinities of 33 to 36 and are on the NBS scale

    KP2 = (exp(-9.039 - 1450 / T)) / fH # convered to SWS pH scale

    KP3 = (exp(4.466 - 7276 / T)) / fH # convered to SWS pH scale

    # Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
    # The Chemical Society (London), Special Publ. 17:751, 1964
    KSi = (0.0000000004 / fH) # convered to SWS pH scale

    return (; KP1=KP1, KP2=KP2, KP3=KP3, KSi=KSi)

end


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
function calc_KP(; T_in, S_in, fH, kwargs...)
    T = T_in + 273.15
    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
    IonS = 19.924 * S_in / (1000 - 1.005 * S_in)

    lnKP1 = (-4576.752 / T + 115.54 - 18.453 * log(T) + (-106.736 / T +
    0.69171) * sqrt(S_in) + (-0.65643 / T - 0.01844) * S_in)
    KP1 = exp(lnKP1)

    lnKP2 = (-8814.715 / T + 172.1033 - 27.927 * log(T) + (-160.34 / T +
    1.3566) * sqrt(S_in) + (0.37335 / T - 0.05778) * S_in)
    KP2 = exp(lnKP2)

    lnKP3 = (-3070.75 / T - 18.126 + (17.27039 / T + 2.81197) * 
    sqrt(S_in) + (-44.99486 / T - 0.09984) * S_in)
    KP3 = exp(lnKP3)

    lnKSi = (-8904.2 / T + 117.4 - 19.334 * log(T) + (-458.79 / T +
    3.5913) * sqrt(IonS) + (188.74 / T - 1.5998) * IonS + (-12.1652 /
    T + 0.07871) * IonS^2)
    KSi = (exp(lnKSi) * (1 - 0.001005109 * S_in))

    return (; KP1=KP1, KP2=KP2, KP3=KP3, KSi=KSi)
end



# Calculate KspA and KspC:
# Mucci, 1983: https://doi.org/10.2475/ajs.283.7.780
function calc_KspA(; T_in, S_in, kwargs...)
    T = T_in + 273.15

    log10KspA = -171.945 - 0.077993 * T + (2903.293 / T) + 71.595 *
    log10(T) + (-0.068393 + 0.0017276 * T + (88.135 / T)) * sqrt(S_in) -
    0.10018 * S_in + 0.0059415 * S_in^(1.5)

    KspA = 10.0^(log10KspA)

    return (; KspA=KspA)
end 


function calc_KspC(; T_in, S_in, kwargs...)
    T = T_in + 273.15

    log10KspC = -171.9065 - 0.077993 * T + (2839.319 / T) + 71.595 *
    log10(T) + (-0.77712 + 0.0028426 * T + (178.34 / T)) * sqrt(S_in) -
    0.07711 * S_in + 0.0041249 * S_in^(1.5)

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
function calc_KH2S(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    lnKH2S = 225.838 - 13275.3 / T - 34.6435 * log(T) + 0.3449 * sqrt(S_in) - 0.0247 * S_in
    KH2S = exp(lnKH2S)
    return KH2S
end


# Calculate KNH3

# NOTES FROM CO2SYS:
# Yao and Millero, Aquatic Geochemistry 1:53-88, 1995   SWS
function Millero_KNH3(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    lnKNH3 = -6285.33 / T + 0.0001635 * T - 0.25444 + (0.46532 - 123.7184 / T) *
    sqrt(S_in) + (-0.01992 + 3.17556 / T) * S_in
    KNH3 = exp(lnKNH3)
    return KNH3
end


# NOTES FROM CO2SYS:
# Clegg Whitfield 1995
# Geochimica et Cosmochimica Acta, Vol. 59, No. 12. pp. 2403-2421
# eq (18)  Total scale   t=[-2 to 40 oC]  S=[0 to 40 ppt]   pK=+-0.00015
function Clegg_KNH3(; T_in, S_in, kwargs...)
    T = T_in + 273.15
    pK = 9.244605 - 2729.33 * (1 / 298.15 - 1 / T)
    pK += (0.04203362 - 11.24742 / T) * S_in^0.25
    pK += (
        -13.6416 + 1.176949 * sqrt(T) - 0.02860785 * T + 545.4834 / T
    ) * sqrt(S_in)
    pK += (
        -0.1462507 + 0.0090226468 * sqrt(T) - 0.0001471361 * T+ 10.5425 / T
    ) * S_in^1.5
    pK += (
        0.004669309 - 0.0001691742 * sqrt(T) - 0.5677934 / T
    ) * S_in^2
    pK += (-2.354039e-05 + 0.009698623 / T) * S_in^2.5
    KNH3 = 10.0^-pK
    # 3. Convert from mol/kg-H2O to mol/kg-SW
    KNH3 = KNH3 * (1 - 0.001005109 * Sal)

    return KNH3
end


# Calculate Ca

# NOTES FROM CO2SYS
# *** CalculateCaforGEOSECS:
# Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
# (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
# Culkin gives Ca = (.0213/40.078)*(Sal/1.80655) in mol/kg-SW
# which corresponds to Ca = .01030*Sal/35.
function Culkin_Ca(; S_in, kwargs...)
    return (0.01026 * S_in / 35) # in mol/kg-SW
end


function RT67_Ca(; S_in, kwargs...)
    return (0.02128 / 40.087 * S_in / 1.80655) # in mol/kg-SW
end


# Helper function to handle when user passes K_method="default"
function which_K(; K_method="default", T_in, S_in,Ca=0.0102821, 
    Mg=0.0528171, kwargs...)
    if Ca != 0.0102821 || Mg != 0.0528171
        return "MyAMI"
    elseif T_in > 35
        return "Millero 2006"
    elseif T_in < 2 && 34.0 < S_in < 37.0
        return "Millero 2002"
    elseif T_in < 0 && 10 < S_in < 50
        return "GP 1989"
    elseif S_in < 1.0
        return "Millero 1979"
    elseif S_in > 50.0
        return "Papadimitriou 2018"
    else 
        return "MyAMI"
    end
end


# Function will be called in calculator.jl in carbon_system, boron_system,
# boron_isotopes, and whole_system to calculate the constants used in ensuing
# calculations.
# Possible options for K_method are "Roy 1993", "GP 1989", "Hansson 1973",
# "DM 1987", "HM 1973", "Mehrbach 1973 A", "Merhbach 1973 B", "Millero 1979",
# "CW 2003", "Lueker 2000", "MPM 2002", "Millero 2002", "Millero 2006",
# "Millero 2010", "Waters 2014", "SB 2020", "Papadimitriou 2018", "Sulpis 2020"
# "MyAMI", or "default".
# If left as "default", user can specify K_mode as "dynamic" or "static". If
# left as default, helper function "which_K" will be called to assess which
# K_method is most appropriate for calculations based on input temperature and
# salinity values. If K_mode="static", then one K_method will be used for the
# entire set of inputs based on average temperature and salinity values.
# Alternatively, if K_mode="dynamic", then which_K will determine the best
# K_method for each sample based on individual salinity and temperature values.
#
# Possible options for KSO4_method are "Dickson", "Khoo", "WM13" or "default".
# Default will be calculated as "Dickson" (reccomended by CO2SYS).
# Possible options for BT_method are "Uppstrom", "Lee", "KSKS18" or "default".
# Default will be calculated as "Uppstrom" (reccomended by CO2SYS).
# Possible options for KF_method are "Dickson", "Perez" or "default".
# Default will be calculated as "Dickson" (reccomended by CO2SYS).
# Possible options for KNH3_method are "Millero", "Clegg" or "default".
# Default will be calculated as "Millero".
# Possible options for Ca_method are "Culkin", "RT67" or "default".
# Default will be calculated as "Culkin".
function K_calculator(; T_in, S_in, P_in=0.0, ST=nothing, FT=nothing,
    BT=nothing, K_method="default", KSO4_method="default", BT_method="default",
    KF_method="default", KNH3_method="default", Ca_method="default",
    K_mode="static", legacy_GEOSECS=false, kwargs...)

    if T_in isa AbstractArray
        
        # If static, determine ONE method for the whole array before calculating
        if K_mode == "static" && K_method == "default"
            T_rep = mean(T_in)
            S_rep = mean(S_in)
            chosen_method = which_K(K_method=K_method, T_in=T_rep, S_in=S_rep; kwargs...)

            raw_results = ((t, s, p) -> K_calculator(
                T_in=t, S_in=s, P_in=p, K_method=chosen_method, K_mode=K_mode; kwargs...
            )).(T_in, S_in, P_in)
            
        else # if dynamic
            raw_results = ((t, s, p) -> K_calculator(
                T_in=t, S_in=s, P_in=p, K_method=K_method, K_mode=K_mode; kwargs...
            )).(T_in, S_in, P_in)
        end

        ##################################
        packed_Ks = (
            K1 = [r.Ks.K1 for r in raw_results],
            K2 = [r.Ks.K2 for r in raw_results],
            K0 = [r.Ks.K0 for r in raw_results],
            KB = [r.Ks.KB for r in raw_results],
            KW = [r.Ks.KW for r in raw_results],
            KS = [r.Ks.KS for r in raw_results],
            KF = [r.Ks.KF for r in raw_results],
            KspA = [r.Ks.KspA for r in raw_results],
            KspC = [r.Ks.KspC for r in raw_results],
            KP1 = [r.Ks.KP1 for r in raw_results],
            KP2 = [r.Ks.KP2 for r in raw_results],
            KP3 = [r.Ks.KP3 for r in raw_results],
            KSi = [r.Ks.KSi for r in raw_results],
            KH2S = [r.Ks.KH2S for r in raw_results],
            KNH3 = [r.Ks.KNH3 for r in raw_results],
        )

        return (
            Ks = packed_Ks,
            ST = [r.ST for r in raw_results],
            FT = [r.FT for r in raw_results],
            BT = [r.BT for r in raw_results],
            Ca = [r.Ca for r in raw_results],
            fH = [r.fH for r in raw_results],
            method = [r.method for r in raw_results] # Changed 'methods' to 'method' for consistency
        )
        #######################################
    end


    if K_method =="default"
        K_method = which_K(K_method=K_method, T_in=T_in, S_in=S_in; kwargs...)
    end

    # --- Boron, Sulfate, Fluoride and Calcium Totals ---
    if BT_method == "Lee"
        final_BT = isnothing(BT) ? Lee_BT(; S_in).BT : BT
    elseif BT_method == "KSK18"
        final_BT = isnothing(BT) ? KSK18_BT(; S_in).BT : BT
    else
        final_BT = isnothing(BT) ? Uppstrom_BT(; S_in).BT : BT
    end

    Ca_in = get(kwargs, :Ca, nothing)

    if Ca_method == ""
        final_Ca = isnothing(Ca_in) ? RT67_Ca(; S_in) : Ca_in
    else
        final_Ca = isnothing(Ca_in) ? Culkin_Ca(; S_in) : Ca_in
    end
    
    final_ST = isnothing(ST) ? calc_ST(; S_in).ST : ST
    final_FT = isnothing(FT) ? calc_FT(; S_in).FT : FT

    # --- Standard KSO4 and KF ---
    KSO4 = if KSO4_method == "Khoo"
        res = Khoo_KS(; T_in, S_in); res isa Number ? res : res.KS
    elseif KSO4_method == "WM13"
        res = WM13_KS(; T_in, S_in); res isa Number ? res : res.KS
    else # Default is Dickson
        res = Dickson_KS(; T_in, S_in); res isa Number ? res : res.KS
    end

    kF = if KF_method == "Perez"
        res = Perez_KF(; T_in, S_in); res isa Number ? res : res.KF
    else
        res = Dickson_KF(; T_in, S_in); res isa Number ? res : res.KF
    end
    kNH3   = KNH3_method == "Clegg" ? Clegg_KNH3(; T_in, S_in) : Millero_KNH3(; T_in, S_in)
    
    # --- Other Baseline Constants ---
    k0   = calc_K0(; T_in, S_in).K0
    fH_val   = calc_fH(; T_in, S_in).fH
    kB   = calc_KB(; T_in, S_in, ST=final_ST, FT=final_FT, KS=KSO4, KF=kF).KB
    kW   = calc_KW(; T_in, S_in).KW
    kPs  = calc_KP(; T_in, S_in, fH=fH_val)
    kH2S = calc_KH2S(; T_in, S_in)
    
    # Initialize K1, K2, KspA, KspC so they exist in memory
    K1 = K2 = KspA = KspC = 1e-10

    # Case 1
    if K_method == "Roy 1993"

        # Calculate K1 & K2
        (; K1, K2) = Roy1993(; T_in, S_in, ST=final_ST, FT=final_FT, KS=KSO4, KF=kF)

    # Case 2
    elseif K_method == "GP 1989"

        # Calculate K1 & K2
        (; K1, K2) = GP1989(; T_in, S_in)

    # Case 3
    elseif K_method == "Hansson 1973"

        # Calculate K1 & K2
        (; K1, K2) = Hansson1973(; T_in, S_in)

    # Case 4
    elseif K_method == "DM 1987"

        # Calculate K1 & K2
        (; K1, K2) = DM1987(; T_in, S_in)

    # Case 5
    elseif K_method == "HM 1973"

        # Calculate K1 & K2
        (; K1, K2) = HM1973(; T_in, S_in)


    # Case 6
    elseif K_method == "Mehrbach 1973 A"
        final_BT = isnothing(BT) ? case67_BT(; S_in).BT : BT
        
        if legacy_GEOSECS
            # Only override these if user explicitly asks for legacy GEOSECS behavior
            final_BT = isnothing(BT) ? case67_BT(; S_in, fH=fH_val).BT : BT
            kB = case67_KB(; T_in, S_in, fH=fH_val).KB
            kW = case6_KW(; T_in, S_in).KW
            kPs = case68_KP(; T_in, S_in, fH=fH_val)
        end
        # Always calculate K1 & K2
        (; K1, K2) = Mehrbach1973(; T_in, S_in, fH=fH_val)

    # Case 7
    elseif K_method == "Mehrbach 1973 B"
        # Calculate fH first so it can be called later
        fH_val = case7_fH(; T_in, S_in).fH
        
        if legacy_GEOSECS
            # Only override these if user explicitly asks for legacy GEOSECS behavior
            final_BT = isnothing(BT) ? case67_BT(; S_in, fH=fH_val).BT : BT
            kB = case67_KB(; T_in, S_in, fH=fH_val).KB
            kW = case7_KW(; T_in, S_in).KW
            kPs = case7_KP(; T_in, S_in, fH=fH_val)
        end
        # Always calculate K1 & K2
        (; K1, K2) = Mehrbach1973(; T_in, S_in, fH=fH_val)


    # Case 8
    elseif K_method == "Millero 1979"

        # Calculate fH first so it can be called later
        fH_val = case8_fH(; T_in, S_in).fH

        # Start by calculating BT, ST and FT if not already passed by user
        final_BT = isnothing(BT) ? case8_BT(; S_in).BT : BT

        # Calculate other constants
        kB = case8_KB(; T_in, S_in, fH_val).KB
        kW = case8_KW(; T_in, S_in).KW
        kPs = case68_KP(; T_in, S_in, fH_val)

        # Calculate K1 & K2
        (; K1, K2) = Millero1979(; T_in, S_in)


    # Case 9
    elseif K_method == "CW 2003"

        # Calculate K1 & K2
        (; K1, K2) = CW2003(; T_in, S_in, fH=fH_val)
    

    # Case 10
    elseif K_method == "Lueker 2000"

        # Calculate K1 & K2
        (; K1, K2) = Lueker2000(; T_in, S_in, ST=final_ST, FT=final_FT, KS=KSO4, KF=kF)


    # Case 11
    elseif K_method == "MPM 2002"

        # Calculate K1 & K2
        (; K1, K2) = MPM2002(; T_in, S_in)

    # Case 12
    elseif K_method == "Millero 2002"

        # Calculate K1 & K2
        (; K1, K2) = Millero2002(; T_in, S_in)

    # Case 13
    elseif K_method == "Millero 2006"

        # Calculate K1 & K2
        (; K1, K2) = Millero2006(; T_in, S_in)


    # Case 14
    elseif K_method == "Millero 2010"

        # Calculate K1 & K2
        (; K1, K2) = Millero2010(; T_in, S_in)

      
    # Case 15
    elseif K_method == "Waters 2014"

        # Calculate K1 & K2
        (; K1, K2) = Waters2014(; T_in, S_in)
    

    # Case 16
    elseif K_method == "SB 2020"

        # Calculate K1 & K2
        (; K1, K2) = SB2020(; T_in, S_in, ST=final_ST, FT=final_FT, KS=KSO4, KF=kF)


    # Case 17
    elseif K_method == "Papadimitriou 2018"

        # Calculate K1 & K2
        (; K1, K2) = Papadimitriou2018(; T_in, S_in, ST=final_ST, FT=final_FT, KS=KSO4, KF=kF)


    # Case # 18
    elseif K_method == "Sulpis 2020"

        # Calculate K1 & K2
        (; K1, K2) = Sulpis2020(; T_in, S_in, ST=final_ST, FT=final_FT, KS=KSO4, KF=kF)


    # MyAMI case
    elseif K_method == "MyAMI"
        
        # Calculate these values based on MyAMI/CBsys equations
        final_ST = isnothing(ST) ? Helpers.calc_ST(S_in) : ST
        final_FT = isnothing(FT) ? Helpers.calc_FT(S_in) : FT
        final_BT = isnothing(BT) ? Helpers.calc_BT(S_in) : BT
        
        Mg_val = get(kwargs, :Mg, 0.0528171) 
        Ca_val = get(kwargs, :Ca, 0.0102821)

        mode_val = get(kwargs, :MyAMI_mode, "calculate")

        py_ks = kgen.calc_Ks(
            temp_c = np.asarray(T_in),
            sal = np.asarray(S_in), 
            p_bar = np.asarray(P_in), 
            magnesium = np.asarray(Mg_val),
            calcium = np.asarray(Ca_val),
            sulphate = np.asarray(final_ST), 
            fluorine = np.asarray(final_FT), 
            MyAMI_mode = mode_val
        )
        
        # Convert to Julia Dict, then to NamedTuple
        ks_dict = pyconvert(Dict{Symbol, Any}, py_ks)
        Ks = (; ks_dict...)

    else
        throw(ArgumentError("Unknown K_method: $K_method"))

    end


    # Accounting for pressure corrections
    P = P_in # in bars
    RGasConstant = 83.14472 # was 83.14462618 ml bar-1 K-1 mol-1, DOEv2, changed to match CO2SYS
    RT = RGasConstant * (T_in + 273.15)

    # Case 8
    if K_method == "Millero 1979" # Correction method from Millero 1983
        deltaV1 = -30.54 + 0.1849 * T_in - 0.0023366 * T_in^2
        Kappa1 = (-5.74 + 0.093 * T_in - 0.001896 * T_in^2) / 1000
        lnK1_fact = (-deltaV1 + 0.5 * Kappa1 * P) * P / RT

        deltaV2 = -29.81 + 0.115 * T_in - 0.001816 * T_in^2
        Kappa2 = (-5.74 + 0.093 * T_in - 0.001896 * T_in^2) / 1000
        lnK2_fact = (-deltaV2 + 0.5 * Kappa2 * P) * P / RT

        lnKB_fact = 0 # this doesn't matter since TB = 0 for this case

    # Cases 6 & 7 (Legacy GEOSECS Pressure)
    elseif (K_method == "Mehrbach 1973 A" || K_method == "Mehrbach 1973 B") && legacy_GEOSECS
        # GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
        lnK1_fact = (24.2 - 0.085 * T_in) * P / RT
        lnK2_fact = (16.4 - 0.04 * T_in) * P / RT
        lnKB_fact = (27.5 - 0.095 * T_in) * P / RT

    # For all others (besides MyAMI, which handles its own pressure corrections)
    elseif !(K_method == "MyAMI")
        # These are from Millero, 1995.
        # They are the same as Millero, 1979 and Millero, 1992.
        # They are from data of Culberson and Pytkowicz, 1968.
        deltaV1 = -25.5 + 0.1271 * T_in
        Kappa1 = (-3.08 + 0.0877 * T_in) / 1000
        lnK1_fact = (-deltaV1 + 0.5 * Kappa1 * P) * P / RT
        # These are from Millero, 1995.
        # They are the same as Millero, 1979 and Millero, 1992.
        # They are from data of Culberson and Pytkowicz, 1968.
        deltaV2 = -15.82 - 0.0219 * T_in
        Kappa2 = (1.13 - 0.1475 * T_in) / 1000
        lnK2_fact = (-deltaV2 + 0.5 * Kappa2 * P) * P / RT
        # This is from Millero, 1979.
        # It is from data of Culberson and Pytkowicz, 1968.
        deltaVB = -29.48 + 0.1622 * T_in - 0.002608 * T_in^2
        KappaB = -2.84 / 1000
        lnKB_fact = (-deltaVB + 0.5 * KappaB * P) * P / RT

    end 

    # Pressure corrections for KW
    # Case 8 (freshwater)
    if K_method == "Millero 1979"
        deltaVW = -25.6 + 0.2324 * T_in - 0.0036246 * T_in^2
        KappaW = (-7.33 + 0.1368 * T_in - 0.001233 * T_in^2) / 1000
        lnKW_fact = ((-deltaVW + 0.5 * KappaW * P) * P / RT)

    elseif K_method ∉ ["MyAMI"]
        deltaVW = -20.02 + 0.1119 * T_in - 0.001409 * T_in^2
        KappaW = (-5.13 + 0.0794 * T_in) / 1000
        lnKW_fact = (-deltaVW + 0.5 * KappaW * P) * P / RT
    end 
    
    # Corrections for KF, KS, KP, and KSi (same for all methods)
    # This is from Millero, 1995, which is the same as Millero, 1983.
    # It is assumed that KF is on the free pH scale.
    deltaVF = -9.78 - 0.009 * T_in - 0.000942 * T_in^2
    KappaF = (-3.91 + 0.054 * T_in) / 1000
    lnKF_fact = (-deltaVF + 0.5 * KappaF * P) * P / RT
    # This is from Millero, 1995, which is the same as Millero, 1983.
    # It is assumed that KS is on the free pH scale.
    deltaVS = -18.03 + 0.0466 * T_in + 0.000316 * T_in^2
    KappaS = (-4.53 + 0.09 * T_in) / 1000
    lnKS_fact = (-deltaVS + 0.5 * KappaS * P) * P / RT
    # The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
    # same as Millero, 1983.
    deltaVP1 = -14.51 + 0.1211 * T_in - 0.000321 * T_in^2
    KappaP1 = (-2.67 + 0.0427 * T_in) / 1000
    lnKP1_fact = (-deltaVP1 + 0.5 * KappaP1 * P) * P / RT

    deltaVP2 = -23.12 + 0.1758 * T_in - 0.002647 * T_in^2
    KappaP2 = (-5.15 + 0.09 * T_in) / 1000
    lnKP2_fact = (-deltaVP2 + 0.5 * KappaP2 * P) * P / RT

    deltaVP3 = -26.57 + 0.202 * T_in - 0.003042 * T_in^2
    KappaP3 = (-4.08 + 0.0714 * T_in) / 1000
    lnKP3_fact = (-deltaVP3 + 0.5 * KappaP3 * P) * P / RT

    deltaVSi = -29.48 + 0.1622 * T_in - 0.002608 * T_in^2
    KappaSi = -2.84 / 1000
    lnKSi_fact = (-deltaVSi + 0.5 * KappaSi * P) * P / RT

    # Corrections for KspA and KspC from Millero 1995:
    deltaVC = -48.76 + 0.5304 * T_in
    KappaC = -11.76e-3 + 0.3692e-3 * T_in
    lnKspC_fact = ((-deltaVC + 0.5 * KappaC * P) * P) / RT

    deltaVA = -46.0 + 0.5304 * T_in
    KappaA = -11.76e-3 + 0.3692e-3 * T_in
    lnKspA_fact = ((-deltaVA + 0.5 * KappaA * P) * P) / RT


    # Pressure corrections for KNH4 and KH2S
    deltaVH2S = -11.07 - 0.009 * T_in - 0.000942 * T_in^2
    KappaH2S = (-2.89 + 0.054 * T_in) / 1000
    lnKH2S_fact = ((-deltaVH2S + 0.5 * KappaH2S * P) * P) / RT

    deltaVNH3 = -26.43 + 0.0889 * T_in - 0.000905 * T_in^2
    KappaNH3 = (-5.03 + 0.0814 * T_in) / 1000
    lnKNH3_fact = ((-deltaVNH3 + 0.5 * KappaNH3 * P) * P) / RT


    # 1. FIRST, calculate the pressure-corrected KS and KF
    KS_press = KSO4 * exp(lnKS_fact)
    KF_press = kF * exp(lnKF_fact)

    # 2. THEN, calculate the sws_factor using those pressure-corrected values
    sws_factor = SWStoTOT(ST=final_ST, FT=final_FT, KS=KS_press, KF=KF_press)

# 3. Determine native pH scales
    sws_native_methods = [
        "GP 1989", "DM 1987", "Hansson 1973",
        "Mehrbach 1973 A", "Mehrbach 1973 B", 
        "MPM 2002", "Millero 2002", "Millero 2006", "Millero 2010",
        "Waters 2014", "HM 1973", "CW 2003"
    ]
    
    nbs_native_methods = ["Millero 1979"]
    KB_sws_native_methods = []

    K1_total = K1
    K2_total = K2
    KB_total = kB
    kH2S_tot = kH2S
    kNH3_tot = kNH3

    # Calculate scale conversion multipliers
    # If SWS, we just need sws_factor. If NBS, we need fH * sws_factor
    function get_scale_mult(method, sws_list, nbs_list, sws_fac, fH)
        if method in sws_list
            return sws_fac
        elseif method in nbs_list
            return sws_fac / fH
        else
            return 1.0 # Already Total
        end
    end

    k1k2_mult = get_scale_mult(K_method, sws_native_methods, nbs_native_methods, sws_factor, fH_val)
    kb_mult   = get_scale_mult(K_method, KB_sws_native_methods, [], sws_factor, fH_val)
    # If Millero KNH3 is SWS, we multiply by sws_factor to get to Total
    knh3_mult = (KNH3_method == "Millero") ? sws_factor : 1.0
    kh2s_mult = 1.0 # Already Total based on your notes

    # 4. Apply pressure corrections, THEN apply the scale shift
    if K_method == "MyAMI"
        final_Ks = merge(Ks, (
            KH2S = kH2S_tot * exp(lnKH2S_fact),
            KNH3 = kNH3_tot * exp(lnKNH3_fact)
        ))
    else
        final_Ks = (
            K1 = (K1_total * exp(lnK1_fact)) * k1k2_mult,
            K2 = (K2_total * exp(lnK2_fact)) * k1k2_mult,
            K0 = k0,
            KB = (KB_total * exp(lnKB_fact)) * kb_mult, 
            KS = KS_press,
            KF = KF_press,
            KH2S = (kH2S_tot * exp(lnKH2S_fact)) * kh2s_mult,
            KW   = (kW * exp(lnKW_fact)) * sws_factor,
            KP1  = (kPs.KP1 * exp(lnKP1_fact)) * sws_factor,
            KP2  = (kPs.KP2 * exp(lnKP2_fact)) * sws_factor,
            KP3  = (kPs.KP3 * exp(lnKP3_fact)) * sws_factor,
            KSi  = (kPs.KSi * exp(lnKSi_fact)) * sws_factor,
            KNH3 = (kNH3 * exp(lnKNH3_fact)) * knh3_mult,
            KspA = calc_KspA(T_in=T_in, S_in=S_in).KspA * exp(lnKspA_fact),
            KspC = calc_KspC(T_in=T_in, S_in=S_in).KspC * exp(lnKspC_fact)
        )
    end

    return (
        Ks = final_Ks,
        ST = final_ST,
        FT = final_FT,
        BT = final_BT,
        Ca = final_Ca,
        fH = fH_val,
        method = K_method
    )

end

export K_calculator
end # module