module Isotopes

include("boron.jl")
using .Boron: calc_chiB

# α fractionation constant & ϵ
"""
Alpha for B fractionation
Klochko, et al., 2006
"""
function get_alphaB()
    return 1.0272
end

"""
Converts alpha (Klochko) to ϵ (which is in delta-space)
"""
function alphaB_to_ϵ(alphaB)
    return (alphaB - 1) * 1000
end


"""
ϵ for B fractionation from Klochko alpha
"""
function get_ϵ()
    return alpha_to_ϵ(get_alphaB())
end


"""
Convert ϵ to alpha
"""
function ϵ_to_alpha(ϵ)
return (ϵ/1000) + 1
end


"""
Converts fractional abundnace (A11) to δ-notation
SRM_ratio: the 11B/10B of SRM, default is NIST951 (4.04367)
(Branson, 2017)
"""
function A11_to_δ11(A11, SRM_ratio=4.04367)
    return ((A11 / (1 - A11) / SRM_ratio - 1)) * 1000
end


"""
Converts fractional abundance (A11) to isotope ratio (R11)
(Branson, 2017)
"""
function A11_to_R11(A11)
    return A11 / (1 - A11)
end


"""
Converts δ-notation (δ11) to fractional abundance (A11)
SRM_ratio: the 11B/10B of SRM, default is NIST951 (4.04367)
(Branson, 2017)
"""
function δ11_to_A11(δ11, SRM_ratio=4.04367)
    return SRM_ratio * (δ11 / 1000 + 1) / (SRM_ratio * (δ11 / 1000 + 1) + 1)
end


"""
Converts δ-notation (δ11) to isotope ratio (R11)
SRM_ratio: the 11B/10B of SRM, default is NIST951 (4.04367)
(Branson, 2017)
"""
function δ11_to_R11(δ11, SRM_ratio=4.04367)
    return (δ11 / 1000 + 1) * SRM_ratio
end


"""
Converts isotope ratio (R11) to δ-notation (δ11)
SRM_ratio: the 11B/10B of SRM, default is NIST951 (4.04367)
(Branson, 2017)
"""
function R11_to_δ11(R11, SRM_ratio=4.04367)
    return (R11 / SRM_ratio - 1) * 1000
end


"""
Converts isotope ratio (R11) to fractional abundance (A11)
(Branson, 2017)
"""
function R11_to_A11(R11)
    return R11 / (1 + R11)
end


"""
Converts the isotope fractional abundnace of B(OH)₃ to isotope fractionation
of B(OH)₄
(Branson, 2017)
"""
function ABOH3_to_ABOH4(ABOH₃, alphaB)
    return (1 / ((alphaB / ABOH₃) - alphaB + 1))
end


"""
Helper function to determne if AB(OH)₃ or AB(OH)₄ is not provided
"""
function ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)
    if all(isnothing, (ABOH₃, ABOH₄))
        throw(ArgumentError("Either AB(OH)₃ or AB(OH)₄ must be specified."))
    elseif isnothing(ABOH₄)
        ABOH₄ = ABOH3_to_ABOH4(ABOH₃, alphaB)
        return ABOH₄
    end
end


"""
Calculates ABT from pH and ABOH₃ or ABOH₄
(Branson, 2017)
"""
function calc_ABT(; H, Ks, alphaB, ABOH₄=nothing, ABOH₃=nothing)
    if isnothing(ABOH₄)
        ABOH₄ = ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)
    end

    chiB = calc_chiB(H, Ks)

    return (
        ABOH₄
        * (
            -ABOH₄ * alphaB * chiB + ABOH₄ * alphaB + ABOH₄ * chiB - ABOH₄
            + alphaB * chiB - chiB + 1
        )
        / (ABOH₄ * alphaB - ABOH₄ + 1)
    )
end


"""
Calculates H⁺ from isotope fractional abundances of boron species
(Branson, 2017)
"""
function H_from_ABOH3_ABOH4(; Ks, alphaB, ABT, ABOH₄=nothing, ABOH₃=nothing)
    if isnothing(ABOH₄)
        ABOH₄ = ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)
    end

    return (Ks.KB / ((alphaB / (1 - ABOH₄ + alphaB * ABOH₄) - 1)
    / (ABT / ABOH₄ - 1) - 1))
end


"""
Calculates AB(OH)₃ from H⁺ and ABT
(Branson, 2017)
"""
function ABOH3_from_H_ABT(H, ABT, Ks, alphaB)
    chiB = calc_chiB(H, Ks)

    return (
        ABT * alphaB - ABT + alphaB * chiB - chiB
        - sqrt(
            ABT^2 * alphaB^2 - 2 * ABT^2 * alphaB + ABT^2
            - 2 * ABT * alphaB^2 * chiB + 2 * ABT * alphaB
            + 2 * ABT * chiB - 2 * ABT + alphaB^2 * chiB^2
            - 2 * alphaB * chiB^2 + 2 * alphaB * chiB + chiB^2
            - 2 * chiB + 1
        )
        + 1
    ) / (2 * chiB * (alphaB - 1))
    
end


"""
Calculates AB(OH)₄ from H⁺ and ABT
(Branson, 2017)
"""
function ABOH4_from_H_ABT(H, ABT, Ks, alphaB)
    chiB = calc_chiB(H, Ks)

    return -(
        ABT * alphaB - ABT - alphaB * chiB + chiB 
        + sqrt(
            ABT^2 * alphaB^2 - 2 * ABT^2 * alphaB 
            + ABT^2 - 2 * ABT * alphaB^2 * chiB + 2 * ABT * alphaB
            + 2 * ABT * chiB - 2 * ABT + alphaB^2 * chiB^2
            - 2 * alphaB * chiB^2 + 2 * alphaB * chiB + chiB^2 - 2 * chiB
            + 1
        )
        - 1 
    ) / (2 * alphaB * chiB - 2 * alphaB - 2 * chiB + 2)
end


"""
Calculates the fractionation factor (alpha) from isotope 
fractionation abundance of ABT and AB(OH)₃
(Branson, 2017)
"""
function alpha_from_ABT_ABOH3(H, Ks, ABT, ABOH₃)
    return ((1
    / ((H / Ks.KB) * (ABT - ABOH₃) + ABT))
    / (ABOH₃ - 1))
end


"""
Calculates the fractionation factor (alpha) form isotope 
fractionation abundance of ABT and AB(OH)₄
(Branson, 2017)
"""
function alpha_from_ABT_ABOH4(H, Ks, ABT, ABOH₄)
    return ((1 / ABOH₄ - 1)
    / (1 / (ABT - ((ABOH₄ - ABT) / (H / Ks. KB))) - 1))
end


"""
Calculates the stoichiometric equilibrium constant for boron
from the fractional abundance of B(OH)₄
(Branson, 2017)
"""
function calc_KB(H, alphaB, ABT, ABOH₄=nothing, ABOH₃=nothing)
    ABOH₄ = ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)

    return (H
    / ((ABOH₄ - ABT) / (ABT - 1 / ((1 / alphaB) * 
    (1 / ABOH₄ - 1) + 1))))
end


"""
Calculates pH, ABT, ABOH₃ and ABOH₄ when two of the four parameters are provided
(Branson, 2017)
"""
function calc_B_isotopes(; pHtot=nothing, ABT=nothing, ABOH₃=nothing, 
    ABOH₄ =nothing, alphaB=nothing, Ks=nothing, kwargs ...)
    if isnothing(alphaB)
        alphaB = get_alphaB() 
    end
    # If pH is known:
    if !isnothing(pHtot)
        H = 10.0^(-pHtot)
    # Use pH to calcultae ABT:
        if isnothing(ABT)
            ABT = calc_ABT(; H, Ks, alphaB, ABOH₄, ABOH₃)
        end
    # If pH is not known and ABT is, use ABT, ABOH₃, and ABOH₄ to calculate:
    else 
        if !isnothing(ABT)
            H = H_from_ABOH3_ABOH4(; Ks, alphaB, ABT, ABOH₄, ABOH₃)
            pHtot = -log10(H)
        else
            throw(ArgumentError(
                "ABT and one of ABOH₃ or ABOH₄ must be specified if pH is missing."))
        end
    end
    # If ABOH₃ is unknown, calculate from H and ABT
    if isnothing(ABOH₃)
        ABOH₃ = ABOH3_from_H_ABT(H, ABT, Ks, alphaB)
    end 
    #If ABOH₄ is unknown, calculate from H and ABT
    if isnothing(ABOH₄)
        ABOH₄ = ABOH4_from_H_ABT(H, ABT, Ks, alphaB)
    end 
    return (; pHtot=pHtot, ABT=ABT, ABOH₃=ABOH₃, ABOH₄=ABOH₄, H=H)
end


"""
Calculates pH on the total scale from δ values
(Branson, 2017)
"""
function pH_from_δ(δ11BT, δ11B4, ϵ=get_ϵ())
    # Calculates fractionation of species from δ values
    ABOH₄ = δ11_to_A11(δ11B4)
    ABT = δ11_to_A11(δ11BT)
    alphaB = ϵ_to_alpha(ϵ)
    # Calculates pH from ABT, ABOH₃, and ABOH₄
    return -log10(H_from_ABOH3_ABOH4(; Ks, alphaB, ABT, ABOH₄, ABOH₃))
end


"""
Calculates pKB from δ inputs
(Branson, 2017)
"""
function pKB_from_δ(pH, δ11BT, δ11B4, ϵ=get_ϵ())
    # Calculates fractionation of species from δ values
    ABOH₄ = δ11_to_A11(δ11B4)
    ABT = δ11_to_A11(δ11BT)
    H = 10.0^(-pH)
    alphaB = ϵ_to_alpha(ϵ)
    # Calculates KB from ABT and ABOH₄
    return -log10(calc_KB(H, alphaB, ABT, ABOH₄))
end


"""
Calculates isotope ratio of total boron in δ units
(Branson, 2017)
"""
function calc_δ11BT(pH, KB, δ11B4, ϵ=get_ϵ())
    ABOH₄ = δ11_to_A11(δ11B4)
    alphaB = ϵ_to_alpha(ϵ)
    H = 10.0^(-pH)

    return A11_to_δ11(calc_ABT(; H, KB, alphaB, ABOH₄))
    
end


"""
Calculates isotope ratio of B(OH)₄ in δ units
(Branson, 2017)
"""
function calc_δ11B4(pH, KB, δ11BT, ϵ=get_ϵ())
    ABT = δ11_toA11(δ11BT)
    alphaB = ϵ_to_alpha(ϵ)

    return A11_to_δ11(ABOH4_from_H_ABT(10.0^(-pH), Ks, ABT, alphaB))
end


"""
Calculates the fractionation factor (ϵ) of B(OH)₃ and B(OH)₄ in δ units
(Branson, 2017)
"""
function calc_ϵ(pH, KB, δ11BT, δ11B4)
    ABOH₄ = δ11_to_A11(ABOH₄)
    ABT = δ11_to_A11(ABT)
    H = 10.0^(-pH)
    alphaB = alpha_from_ABT_ABOH4(H, KB, ABT, ABOH₄)

    return alpha_to_ϵ(alphaB)
    
end
export calc_B_isotopes, A11_to_δ11
end # module