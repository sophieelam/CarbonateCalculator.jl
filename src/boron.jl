"""
    Boron

Aqueous boron speciation: the partitioning of total boron between B(OH)₃ and B(OH)₄⁻.

`KB` governs `B(OH)₃ + H₂O ⇌ B(OH)₄⁻ + H⁺`, so a single constraint beyond `BT` fixes the
system. [`B_calculator`](@ref) is the entry point; the rest are the individual rearrangements
it dispatches to.

Equations follow Branson (2017), *B systematics in cbsyst*.

!!! warning
    These functions take concentrations in mol/kg and `Ks` as an equilibrium-constant bundle,
    not the units a result is reported in. Passing values off `result.val`, which are in the
    reporting unit (µmol/kg by default), gives a plausible wrong answer.
"""
module Boron

"""
    calc_chiB(H, Ks)

Return the fraction of total boron present as B(OH)₃, at a given [H⁺] in mol/kg.
"""
function calc_chiB(H, Ks)
    return 1 / (1 + Ks.KB / H)
end

"""
    H_from_BT_BOH3(BT, BOH₃, Ks)

Return [H⁺] in mol/kg, from total boron and B(OH)₃.
"""
function H_from_BT_BOH3(BT, BOH₃, Ks)
    return Ks.KB / (BT / BOH₃ - 1)
end


"""
    H_from_BT_BOH4(BT, BOH₄, Ks)

Return [H⁺] in mol/kg, from total boron and B(OH)₄⁻.
"""
function H_from_BT_BOH4(BT, BOH₄, Ks)
    return Ks.KB * (BT / BOH₄ - 1)
end


"""
    BT_from_pH_BOH3(pH, BOH₃, Ks)

Return total boron in mol/kg, from pH on the total scale and B(OH)₃.
"""
function BT_from_pH_BOH3(pH, BOH₃, Ks)
    H = 10.0^(-pH)
    return BOH₃ * (1 + Ks.KB / H)
end


"""
    BT_from_pH_BOH4(pH, BOH₄, Ks)

Return total boron in mol/kg, from pH on the total scale and B(OH)₄⁻.
"""
function BT_from_pH_BOH4(pH, BOH₄, Ks)
    H = 10.0^(-pH)
    return BOH₄ * (1 + H / Ks.KB)
end



"""
    calc_BOH4(BT, H, Ks)

Return B(OH)₄⁻ in mol/kg, from total boron and [H⁺].
"""
function calc_BOH4(BT, H, Ks)
    return BT / (1 + H / Ks.KB)
end


"""
    calc_BOH3(BT, H, Ks)

Return B(OH)₃ in mol/kg, from total boron and [H⁺].
"""
function calc_BOH3(BT, H, Ks)
    return BT / (1 + Ks.KB / H)
end


"""
    B_calculator(; pHtot, BT, BOH₃, BOH₄, Ks)

Solve boron speciation from any two of `pHtot`, `BT`, `BOH₃` and `BOH₄`.

Concentrations are in mol/kg and pH is on the total scale. Returns
`(; pHtot, BT, BOH₃, BOH₄)`, with the two that were not supplied filled in.
"""
function B_calculator(; pHtot=nothing, BT=nothing, BOH₃=nothing, BOH₄=nothing,
    Ks=nothing, kwargs...)

    # If pH and BT are known, convert pH to H⁺:
    if !isnothing(pHtot) && !isnothing(BT)
        H = 10.0^(-pHtot)
    # If BT and BOH₃ are known, calculate H⁺:
    elseif !isnothing(BT) && !isnothing(BOH₃)
        H = H_from_BT_BOH3(BT, BOH₃, Ks)
    # If BT and BOH₄ are known, calculate H⁺:
    elseif !isnothing(BT) && !isnothing(BOH₄)
        H = H_from_BT_BOH4(BT, BOH₄, Ks)
    # If B(OH)₃ and B(OH)₄ are known, calculate BT and H⁺:
    elseif !isnothing(BOH₃) && !isnothing(BOH₄)
        BT = BOH₃ + BOH₄
        H = H_from_BT_BOH3(BT, BOH₃, Ks)
    # If pH and BOH₃ are known, calculate H⁺ and BT:
    elseif !isnothing(pHtot) && !isnothing(BOH₃)
        H = 10.0^(-pHtot)
        BT = BT_from_pH_BOH3(pHtot, BOH₃, Ks)
    # If pH and BOH₄ are known, calculate H⁺ and BT:
    elseif !isnothing(pHtot) && !isnothing(BOH₄)
        H = 10.0^(-pHtot)
        BT = BT_from_pH_BOH4(pHtot, BOH₄, Ks)
    end 
    
    # Above calculations ensure [H⁺] and BT are known, allowing the
    # remaining species to be calculated.

    if isnothing(BOH₃)
        BOH₃ = calc_BOH3(BT, H, Ks)
    end

    if isnothing(BOH₄)
        BOH₄ = calc_BOH4(BT, H, Ks)
    end

    if isnothing(pHtot)
        pHtot = -log10(H)
    end

    return (pHtot=pHtot, BT=BT, BOH₃=BOH₃, BOH₄=BOH₄)

end 
export B_calculator
end # module
