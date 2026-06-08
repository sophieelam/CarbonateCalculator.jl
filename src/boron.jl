module Boron

"""
Calculation for Chi B
(Branson, 2017)
"""
function calc_chiB(H, Ks)
    return 1 / (1 + Ks.KB / H)
end

"""
#1: Calculating [H⁺] from BT and B(OH)₃
Taken from CBsyst (Branson, 2017)
"""
function H_from_BT_BOH3(BT, BOH₃, Ks)
    return Ks.KB / (BT / BOH₃ - 1)
end


"""
#2: Calculating [H⁺] from BT and B(OH)₄
Taken from CBsyst (Branson, 2017)
"""
function H_from_BT_BOH4(BT, BOH₄, Ks)
    return Ks.KB * (BT / BOH₄ - 1)
end


"""
#3: Calculating BT from pH and B(OH)₃
Taken from CBsyst (Branson, 2017)
"""
function BT_from_pH_BOH3(pH, BOH₃, Ks)
    H = 10.0^(-pH)
    return BOH₃ * (1 + Ks.KB / H)
end


"""
#4: Calculating BT from pH and B(OH)₄
Taken from CBsyst (Branson, 2017)
"""
function BT_from_pH_BOH4(pH, BOH₄, Ks)
    H = 10.0^(-pH)
    return BOH₄ * (1 + H / Ks.KB)
end



"""
Calculating B(OH)₄ from BT and H⁺
Taken from CBsyst (Branson, 2017)
"""
function calc_BOH4(BT, H, Ks)
    return BT / (1 + H / Ks.KB)
end


"""
Calculating B(OH)₃ from BT and H⁺
Taken from CBsyst (Branson, 2017)
"""
function calc_BOH3(BT, H, Ks)
    return BT / (1 + Ks.KB / H)
end


"""
Calculates boron system from any two of the following: 
pH, BT, B(OH)₃, B(OH)₄
(Branson, 2017)
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
