"""
    Isotopes

Boron isotopes and the fractionation between the two dissolved species.

¹¹B partitions unevenly between B(OH)₃ and B(OH)₄⁻, by a factor `alphaB`, so the isotopic
composition of either species is a function of pH.

Two notations appear throughout, and they are the same quantity written differently: `A`, the
fractional abundance of ¹¹B, is what the equations use, and `δ`, the per-mil deviation from a
standard, is what people measure. The `A11_to_δ11` family converts between them, along with
the raw isotope ratio `R`.

[`calc_B_isotopes`](@ref) is the entry point. Equations follow Branson (2017),
*B systematics in cbsyst*; `alphaB` is from Klochko et al. (2006).

!!! warning
    These functions take concentrations in mol/kg and `Ks` as an equilibrium-constant bundle,
    not the units a result is reported in. Fractional abundances are dimensionless and δ
    values are in ‰.
"""
module Isotopes

using ..Boron: calc_chiB

# α fractionation constant & ϵ
"""
    get_alphaB()

Return the boron isotope fractionation factor between B(OH)₃ and B(OH)₄⁻, dimensionless.

Klochko et al., 2006, doi:10.1016/j.epsl.2006.05.034.
"""
function get_alphaB()
    return 1.0272
end

"""
    get_δBT()

Return δ¹¹B of modern seawater, in ‰.

Used whenever the total boron isotope composition is not given explicitly.

Foster et al., 2010, doi:10.1029/2010GC003201.
"""
function get_δBT()
    return 39.61
end

"""
    alphaB_to_ϵ(alphaB)

Return the fractionation factor as ϵ in ‰, the form that works in δ-space.
"""
function alphaB_to_ϵ(alphaB)
    return (alphaB - 1) * 1000
end


"""
    get_ϵ()

Return the boron isotope fractionation factor as ϵ in ‰.

[`get_alphaB`](@ref) expressed in δ-space.
"""
function get_ϵ()
    return alphaB_to_ϵ(get_alphaB())
end


"""
    ϵ_to_alpha(ϵ)

Return the fractionation factor as a dimensionless alpha, from ϵ in ‰.
"""
function ϵ_to_alpha(ϵ)
return (ϵ/1000) + 1
end


"""
    A11_to_δ11(A11, SRM_ratio = 4.04367)

Return δ¹¹B in ‰, from the fractional abundance of ¹¹B.

`SRM_ratio` is the ¹¹B/¹⁰B of the standard, NIST951 by default.
"""
function A11_to_δ11(A11, SRM_ratio=4.04367)
    return ((A11 / (1 - A11) / SRM_ratio - 1)) * 1000
end


"""
    A11_to_R11(A11)

Return the ¹¹B/¹⁰B ratio, from the fractional abundance of ¹¹B.
"""
function A11_to_R11(A11)
    return A11 / (1 - A11)
end


"""
    δ11_to_A11(δ11, SRM_ratio = 4.04367)

Return the fractional abundance of ¹¹B, from δ¹¹B in ‰.

`SRM_ratio` is the ¹¹B/¹⁰B of the standard, NIST951 by default.
"""
function δ11_to_A11(δ11, SRM_ratio=4.04367)
    return SRM_ratio * (δ11 / 1000 + 1) / (SRM_ratio * (δ11 / 1000 + 1) + 1)
end


"""
    δ11_to_R11(δ11, SRM_ratio = 4.04367)

Return the ¹¹B/¹⁰B ratio, from δ¹¹B in ‰.

`SRM_ratio` is the ¹¹B/¹⁰B of the standard, NIST951 by default.
"""
function δ11_to_R11(δ11, SRM_ratio=4.04367)
    return (δ11 / 1000 + 1) * SRM_ratio
end


"""
    R11_to_δ11(R11, SRM_ratio = 4.04367)

Return δ¹¹B in ‰, from the ¹¹B/¹⁰B ratio.

`SRM_ratio` is the ¹¹B/¹⁰B of the standard, NIST951 by default.
"""
function R11_to_δ11(R11, SRM_ratio=4.04367)
    return (R11 / SRM_ratio - 1) * 1000
end


"""
    R11_to_A11(R11)

Return the fractional abundance of ¹¹B, from the ¹¹B/¹⁰B ratio.
"""
function R11_to_A11(R11)
    return R11 / (1 + R11)
end


"""
    ABOH3_to_ABOH4(ABOH₃, alphaB)

Return the fractional abundance of ¹¹B in B(OH)₄⁻, from that in B(OH)₃.
"""
function ABOH3_to_ABOH4(ABOH₃, alphaB)
    return (1 / ((alphaB / ABOH₃) - alphaB + 1))
end


"""
    ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)

Return AB(OH)₄ given either species, converting from AB(OH)₃ when AB(OH)₄ is absent.

Raises when neither is given, so a caller may assign the result unconditionally.
"""
function ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)
    if all(isnothing, (ABOH₃, ABOH₄))
        throw(ArgumentError("Either AB(OH)₃ or AB(OH)₄ must be specified."))
    elseif isnothing(ABOH₄)
        return ABOH3_to_ABOH4(ABOH₃, alphaB)
    end
    return ABOH₄
end


"""
    calc_ABT(; H, Ks, alphaB, ABOH₄, ABOH₃)

Return the fractional abundance of ¹¹B in total boron, from [H⁺] and one of the two species.
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
    H_from_ABOH3_ABOH4(; Ks, alphaB, ABT, ABOH₄, ABOH₃)

Return [H⁺] in mol/kg, from total boron's isotopic composition and one of the two species.

This is the step that makes δ¹¹B a pH proxy.
"""
function H_from_ABOH3_ABOH4(; Ks, alphaB, ABT, ABOH₄=nothing, ABOH₃=nothing)
    if isnothing(ABOH₄)
        ABOH₄ = ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)
    end

    return (Ks.KB / ((alphaB / (1 - ABOH₄ + alphaB * ABOH₄) - 1)
    / (ABT / ABOH₄ - 1) - 1))
end


"""
    ABOH3_from_H_ABT(H, ABT, Ks, alphaB)

Return the fractional abundance of ¹¹B in B(OH)₃, from [H⁺] and total boron.
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
    ABOH4_from_H_ABT(H, ABT, Ks, alphaB)

Return the fractional abundance of ¹¹B in B(OH)₄⁻, from [H⁺] and total boron.
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
    alpha_from_ABT_ABOH3(H, Ks, ABT, ABOH₃)

Return the fractionation factor, from [H⁺] and the composition of total boron and B(OH)₃.
"""
function alpha_from_ABT_ABOH3(H, Ks, ABT, ABOH₃)
    return ((1
    / ((H / Ks.KB) * (ABT - ABOH₃) + ABT))
    / (ABOH₃ - 1))
end


"""
    alpha_from_ABT_ABOH4(H, Ks, ABT, ABOH₄)

Return the fractionation factor, from [H⁺] and the composition of total boron and B(OH)₄⁻.
"""
function alpha_from_ABT_ABOH4(H, Ks, ABT, ABOH₄)
    return ((1 / ABOH₄ - 1)
    / (1 / (ABT - ((ABOH₄ - ABT) / (H / Ks. KB))) - 1))
end


"""
    calc_KB(H, alphaB, ABT, ABOH₄, ABOH₃)

Return the stoichiometric equilibrium constant for boron, from the isotopic composition.

Takes either species; give `ABOH₃` when `ABOH₄` is unknown.
"""
function calc_KB(H, alphaB, ABT, ABOH₄=nothing, ABOH₃=nothing)
    ABOH₄ = ABOH3_or_ABOH4(ABOH₃, ABOH₄, alphaB)

    return (H
    / ((ABOH₄ - ABT) / (ABT - 1 / ((1 / alphaB) * 
    (1 / ABOH₄ - 1) + 1))))
end


"""
    calc_B_isotopes(; pHtot, ABT, ABOH₃, ABOH₄, alphaB, Ks)

Solve the boron isotope system from `pHtot`, or from `ABT` and one of the two species.

Isotopic compositions are fractional abundances of ¹¹B, and pH is on the total scale.
`alphaB` falls back to [`get_alphaB`](@ref). Returns `(; pHtot, ABT, ABOH₃, ABOH₄, H)`.
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
    pH_from_δ(KB, δ11BT, δ11B4, ϵ = get_ϵ())

Return pH on the total scale, from the δ¹¹B of total boron and of B(OH)₄⁻, both in ‰.

This is the boron isotope pH proxy. `KB` is the stoichiometric constant for boron and `ϵ` the
fractionation in ‰.

Inverts [`calc_δ11B4`](@ref).
"""
function pH_from_δ(KB, δ11BT, δ11B4, ϵ=get_ϵ())
    ABOH₄ = δ11_to_A11(δ11B4)
    ABT = δ11_to_A11(δ11BT)
    alphaB = ϵ_to_alpha(ϵ)

    return -log10(H_from_ABOH3_ABOH4(; Ks = (; KB), alphaB, ABT, ABOH₄))
end


"""
    pKB_from_δ(pH, δ11BT, δ11B4, ϵ = get_ϵ())

Return pKB, from pH on the total scale and the δ¹¹B of total boron and of B(OH)₄⁻, in ‰.

The inverse of [`pH_from_δ`](@ref): given an independently known pH, what `KB` the isotopes
imply.
"""
function pKB_from_δ(pH, δ11BT, δ11B4, ϵ=get_ϵ())
    ABOH₄ = δ11_to_A11(δ11B4)
    ABT = δ11_to_A11(δ11BT)
    H = 10.0^(-pH)
    alphaB = ϵ_to_alpha(ϵ)

    return -log10(calc_KB(H, alphaB, ABT, ABOH₄))
end


"""
    calc_δ11BT(pH, KB, δ11B4, ϵ = get_ϵ())

Return δ¹¹B of total boron in ‰, from pH on the total scale and the δ¹¹B of B(OH)₄⁻.

Inverts [`calc_δ11B4`](@ref).
"""
function calc_δ11BT(pH, KB, δ11B4, ϵ=get_ϵ())
    ABOH₄ = δ11_to_A11(δ11B4)
    alphaB = ϵ_to_alpha(ϵ)
    H = 10.0^(-pH)

    return A11_to_δ11(calc_ABT(; H, Ks = (; KB), alphaB, ABOH₄))
end


"""
    calc_δ11B4(pH, KB, δ11BT, ϵ = get_ϵ())

Return δ¹¹B of B(OH)₄⁻ in ‰, from pH on the total scale and the δ¹¹B of total boron.

The inverse of [`calc_δ11BT`](@ref), and the quantity a carbonate archive records.
"""
function calc_δ11B4(pH, KB, δ11BT, ϵ=get_ϵ())
    ABT = δ11_to_A11(δ11BT)
    alphaB = ϵ_to_alpha(ϵ)

    return A11_to_δ11(ABOH4_from_H_ABT(10.0^(-pH), ABT, (; KB), alphaB))
end


"""
    calc_ϵ(pH, KB, δ11BT, δ11B4)

Return the fractionation factor ϵ in ‰, from pH and the δ¹¹B of total boron and of B(OH)₄⁻.

Solves for the fractionation a measurement implies, rather than assuming
[`get_ϵ`](@ref) — the calibration question behind the proxy.
"""
function calc_ϵ(pH, KB, δ11BT, δ11B4)
    ABOH₄ = δ11_to_A11(δ11B4)
    ABT = δ11_to_A11(δ11BT)
    H = 10.0^(-pH)
    alphaB = alpha_from_ABT_ABOH4(H, (; KB), ABT, ABOH₄)

    return alphaB_to_ϵ(alphaB)
end
export calc_B_isotopes, A11_to_δ11
end # module