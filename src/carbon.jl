module Carbon
## Remember to run packages!
using Roots
using ForwardDiff
include("helpers.jl")
using .Helpers

export C_calculator, calc_revelle_factor, calc_buffer_capacity, fCO₂_to_CO₂,
CO₂_to_fCO₂, fCO₂_to_pCO₂, pCO₂_to_fCO₂

## Calculation #1
"""
#1: Calculating DIC from CO₂ and pH
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""

function DIC_from_CO₂_pH(CO₂, pH, Ks)
    H = 10.0^(-pH)
    return CO₂ * (1.0 + Ks.K1/H + Ks.K1*Ks.K2/H^2)
end


## Calculation #2 
"""
#2: Calculating H⁺ from CO₂ and HCO₃⁻
Zeebe & Wolf-Gladrow, 2001, Appendix B
Solved using autograd
"""
# function solve_H_from_CO₂_HCO₃(H, CO₂, HCO₃, Ks)
#     LH = CO₂ * (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2)
#     RH = HCO₃ * (H^2 + H^3 / Ks.K1 + Ks.K2 * H)
#     return LH - RH
# end


# function H_from_CO₂_HCO₃(CO₂, HCO₃, Ks)
#     f(H) = solve_H_from_CO₂_HCO₃(H, CO₂, HCO₃, Ks)
#     df(H) = ForwardDiff.derivative(f, H)
    
#     initial_guess = 1e-8 + + zero(CO₂) + zero(HCO₃) # Standard starting guess (~pH 8)
#     return find_zero((f, df), initial_guess, Roots.Newton())
# end

function H_from_CO₂_HCO₃(CO₂, HCO₃, Ks)
    return (Ks.K1 * CO₂) / HCO₃
end

## Calculation #3 
"""
#3: Calculating H⁺ from CO₂ and CO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
Solved using autograd
"""
# function solve_H_from_CO₂_CO₃(H, CO₂, CO₃, Ks)
#     LH = CO₂ * (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2)
#     RH = CO₃ * (H^2 + H^3 / Ks.K2 + H^4 / (Ks.K1 * Ks.K2))
#     return LH - RH
# end

# function H_from_CO₂_CO₃(CO₂, CO₃, Ks)
#     f(H) = solve_H_from_CO₂_CO₃(H, CO₂, CO₃, Ks)
#     df(H) = ForwardDiff.derivative(f, H)
    
#     initial_guess = 1e-8 + zero(CO₂) + zero(CO₃)
#     return find_zero((f, df), initial_guess, Roots.Newton())
# end

function H_from_CO₂_CO₃(CO₂, CO₃, Ks)
    # Using abs() just as a safety net against tiny floating-point noise around zero
    return sqrt(abs((Ks.K1 * Ks.K2 * CO₂) / CO₃))
end


## Calculation #4 
"""
#4: Calculating pH from CO₂ and TA
Taken from MatLab CO2SYS (which originally used a Newton-Raphson method) and
adapted to be solved more efficiently with the Julia ForwardDiff auto-grad capabilites.
"""
function solve_pH_from_CO₂_TA(pH, CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

    H = 10.0^(-pH)
    fCO₂ = CO₂ / Ks.K0
    HCO₃ = Ks.K0 * Ks.K1 * fCO₂ / H
    CO₃ = Ks.K0 * Ks.K1 * Ks.K2 * fCO₂ / H^2
    CAlk = HCO₃ + 2 * CO₃
    BAlk = BT * Ks.KB / (Ks.KB + H)
    OH = Ks.KW / H
    PhosNum = Ks.KP1 * Ks.KP2 * H + 2 * Ks.KP1 * Ks.KP2 * Ks.KP3 - H^3
    PhosDenom = H^3 + Ks.KP1 * H^2 + Ks.KP1 * Ks.KP2 * H + Ks.KP1 * Ks.KP2 * Ks.KP3
    PAlk = PT * PhosNum / PhosDenom
    SiAlk = SiT * Ks.KSi / (Ks.KSi + H)
    Alk_H2S = H2ST * Ks.KH2S / (Ks.KH2S + H)
    Alk_NH3 = NH4T * Ks.KNH3 / (Ks.KNH3 + H)
    Hfree = H / (1 + ST / Ks.KS)
    HSO₄ = ST / (1 + Ks.KS / Hfree)
    HF = FT / (1 + Ks.KF / Hfree)

    return TA - CAlk - BAlk - OH - PAlk - SiAlk - Alk_H2S - Alk_NH3 + Hfree + HSO₄ + HF
end

function pH_from_CO₂_TA(CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # 1. Create temporary function
    f(pH) = solve_pH_from_CO₂_TA(pH, CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    
    # 2. Use ForwardDiff for exact derivative
    df(pH) = ForwardDiff.derivative(f, pH)
    
    # 3. Solve using Newton's method
    initial_guess = 8.0 + zero(CO₂) + zero(TA)
    return find_zero((f, df), initial_guess, Roots.Newton())
end


# function H_from_HCO₃_CO₃(HCO₃, CO₃, Ks)
#     f(H) = solve_H_from_HCO₃_CO₃(H, HCO₃, CO₃, Ks)
#     df(H) = ForwardDiff.derivative(f, H)
    
#     initial_guess = 1e-8 + zero(HCO₃) + zero(CO₃)
#     return find_zero((f, df), initial_guess, Roots.Newton())
# end

# function H_from_HCO₃_CO₃(HCO₃, CO₃, Ks)
#     return (Ks.K2 * HCO₃) / CO₃
# end


## Calculation #5 
"""
#5: Calculating H⁺ from CO₂ and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B

"""
function solve_H_from_CO₂_DIC(H, CO₂, DIC, Ks)
    LH = DIC * H^2
    RH = CO₂ * (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2)
    return LH - RH
end


function H_from_CO₂_DIC(CO₂, DIC, Ks)
    f(H) = solve_H_from_CO₂_DIC(H, CO₂, DIC, Ks)
    return find_zero(f, (1e-14, 1))
end


## Calculation #6
"""
#6: Calculating DIC from pH and HCO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function DIC_from_pH_HCO₃(pH, HCO₃, Ks)
    H = 10.0^(-pH)
    return HCO₃ * (1.0 + H / Ks.K1 + Ks.K2 / H)
end


## Calculation #7
"""
#7: Calculating DIC from pH and CO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function DIC_from_pH_CO₃(pH, CO₃, Ks)
    H = 10.0^(-pH)
    return CO₃ * (1 + H / Ks.K2 + H^2 / (Ks.K1 * Ks.K2))
end


## Calculation #8
"""
#8: Calculating DIC from pH and TA
Taken from MatLab CO2SYS
"""
function DIC_from_pH_TA(pH, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
   H = 10 ^(-pH)
   BAlk = BT * Ks.KB / (Ks.KB + H)
    OH = Ks.KW / H
    PhosNum = Ks.KP1 * Ks.KP2 * H + 2 * Ks.KP1 * Ks.KP2 * Ks.KP3 - H^3
    PhosDenom = H^3 + Ks.KP1 * H^2 + Ks.KP1 * Ks.KP2 * H + Ks.KP1 * Ks.KP2 * Ks.KP3
    PAlk = PT * PhosNum / PhosDenom
    SiAlk = SiT * Ks.KSi / (Ks.KSi + H)
    Alk_H2S = H2ST * Ks.KH2S / (Ks.KH2S + H)
    Alk_NH3 = NH4T * Ks.KNH3 / (Ks.KNH3 + H)
    Hfree = H / (1 + ST / Ks.KS)
    HSO₄ = ST / (1 + Ks.KS / Hfree)
    HF = FT / (1 + Ks.KF / Hfree)
    CAlk = TA - BAlk - OH - PAlk - SiAlk - Alk_H2S - Alk_NH3 + Hfree + HSO₄ + HF
    return CAlk * (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2) / (Ks.K1 * (H + 2.0 * Ks.K2))
end

## Calculation #9
"""
#9: Calculating CO₂ from pH and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function CO₂_from_pH_DIC(pH, DIC, Ks)
    H = 10.0^(-pH)
    return DIC / (1.0 + Ks.K1 / H + Ks.K1 * Ks.K2 / H^2)
end


## Calculation #10
"""
#10: Calculating H⁺ from HCO₃ and CO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function solve_H_from_HCO₃_CO₃(H, HCO₃, CO₃, Ks)
    LH = HCO₃ * (H + H^2 / Ks.K1 + Ks.K2)
    RH = CO₃ * (H + H ^2 / Ks.K2 + H^3 / (Ks.K1 * Ks.K2))
    return LH - RH
end


function H_from_HCO₃_CO₃(HCO₃, CO₃, Ks)
    f(H) = solve_H_from_HCO₃_CO₃(H, HCO₃, CO₃, Ks)
    return find_zero(f, (1e-14, 1))
end 



## Calculation #11
"""
#11: Calculating H⁺ from HCO₃ and TA
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function solve_H_from_HCO₃_TA(H, HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    temp_DIC = HCO₃ * (H / Ks.K1 + 1.0 + Ks.K2 / H)
    (calc_TA_val, _, _, _, _, _, _, _, _) = calc_TA(H, temp_DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks; mode="multi")
    return calc_TA_val - TA
end

function H_from_HCO₃_TA(HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    f(pH) = solve_H_from_HCO₃_TA(10.0^(-pH), HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    df(pH) = ForwardDiff.derivative(f, pH)
    initial_guess = 8.0 + zero(HCO₃) + zero(TA)
    sol_pH = find_zero((f, df), initial_guess, Roots.Newton())
    return 10.0^(-sol_pH)
end

# function solve_H_from_HCO₃_TA(H, HCO₃, TA, BT, Ks)
#     LH = TA * (Ks.KB + H) * (H^3 + Ks.K1 * H^2 + Ks.K1 * Ks.K2 * H)
#     RH = (
#         HCO₃ * (H + H^2 / Ks.K1 + Ks.K2) 
#         * ((Ks.KB + 2 * Ks.K2) * Ks.K1 * H + 2 * Ks.KB * Ks.K1 * Ks.K2 + Ks.K1 * H^2)
#         + ((H^2 + Ks.K1 * H + Ks.K1 * Ks.K2)
#         *(Ks.KB * BT * H + Ks.KW * Ks.KB + Ks.KW * H - Ks.KB * H^2 - H^3))
#     )
#     return LH - RH
# end

# function H_from_HCO₃_TA(HCO₃, TA, BT, Ks)
#     f(H) = solve_H_from_HCO₃_TA(H, HCO₃, TA, BT, Ks)
#     return find_zero(f, (1e-14, 1))
# end 


## Calculation #12
"""
#12: Calculating pH from HCO₃ and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B
Note: instead of using "find_zero", this calculation uses a basic quadratic 
approach to save computational time
"""
function pH_from_HCO₃_DIC(HCO₃, DIC, Ks)
    # Rearranging the equation into aH^2 + bH + c = 0
    a = HCO₃ / Ks.K1
    b = HCO₃ - DIC
    c = HCO₃ * Ks.K2
    
    discriminant = b^2 - 4 * a * c
    
    if discriminant < 0
        return NaN
    end
    
    # Calculate both mathematical roots
    H_1 = (-b + sqrt(discriminant)) / (2 * a)
    H_2 = (-b - sqrt(discriminant)) / (2 * a)
    
    larger_H = maximum([H_1, H_2])
    smaller_H = minimum([H_1, H_2])
    
    # Returns the smaller H⁺ concentration of two roots found (higher pH)
    return -log10(smaller_H) # This doesn't match CBsyst, but matches test suite
end


## Calculation #13
"""
#13: Calculating H⁺ from CO₃ and TA
Zeebe & Wolf-Gladrow, 2001, Appendix B
Uses Roots.Brent() (same as CBsys) to circumnavigate bracketing issues with
root finding. However, this only works for pH values 5 < pH < 10.
"""
function solve_H_from_CO₃_TA(H, CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    temp_DIC = CO₃ * (H^2 / (Ks.K1 * Ks.K2) + H / Ks.K2 + 1)
    (calc_TA_val, _, _, _, _, _, _, _, _) = calc_TA(H, temp_DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks; mode="multi")
    return calc_TA_val - TA
end

function H_from_CO₃_TA(CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # 1. Create stability wrapper
    f(pH) = solve_H_from_CO₃_TA(10.0^(-pH), CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    
    # 2. Get derivative
    df(pH) = ForwardDiff.derivative(f, pH)
    
    # 3. Solve
    initial_guess = 8.0 + zero(CO₃) + zero(TA)
    sol_pH = find_zero((f, df), initial_guess, Roots.Newton())
    
    # 4. Return H⁺
    return 10.0^(-sol_pH)
end

# function solve_H_from_CO₃_TA(H, CO₃, TA, BT, Ks)
#     LH = TA * (Ks.KB + H) * (H^3 + Ks.K1 * H^2 + Ks.K1 * Ks.K2 * H)
#     RH = (
#         CO₃ * (H + H^2 / Ks.K2 + H^3 / (Ks.K1 * Ks.K2))
#         * (Ks.K1 * H^2 + Ks.K1 * H * (Ks.KB + 2 * Ks.K2) + 2 * Ks.KB * Ks.K1 * Ks.K2)
#     ) + (
#         (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2)
#         * (Ks.KB * BT * H + Ks.KW * Ks.KB + Ks.KW * H - Ks.KB * H^2 - H^3) 
#     )
    
#     return LH - RH
# end

# function H_from_CO₃_TA(CO₃, TA, BT, Ks)
#     f(H) = solve_H_from_CO₃_TA(H, CO₃, TA, BT, Ks)
    
#     # Expanded bracket for GLODAP, and explicitly matching Python's algorithm
#     return find_zero(f, (1e-10, 1e-5), Roots.Brent()) 
# end


## Calculation #14
"""
#14: Calculating H⁺ from CO₃ and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function solve_H_from_CO₃_DIC(H, CO₃, DIC, Ks)
    LH = CO₃ * (1 + H / Ks.K2 + H^2 / (Ks.K1 * Ks.K2))
    RH = DIC
    return LH - RH
end

function H_from_CO₃_DIC(CO₃, DIC, Ks)
    f(H) = solve_H_from_CO₃_DIC(H, CO₃, DIC, Ks)
    return find_zero(f, (1e-14, 1))
end


## Calculation #15 
"""
#15: Calculating pH from TA and DIC
Taken from MatLab CO2SYS (which originally used a Newton-Raphson method) and
adapted to be solved more efficiently with Julia ForwardDiff autograd capabilites.
"""
function solve_pH_from_TA_DIC(pH, TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

    H = 10.0^(-pH)
    Denom = H^2 + Ks.K1 * H + Ks.K1 * Ks.K2
    CAlk = DIC * Ks.K1 * (H + 2 * Ks.K2) / Denom
    BAlk = BT * Ks.KB / (Ks.KB + H)
    OH = Ks.KW / H
    PhosNum = Ks.KP1 * Ks.KP2 * H + 2 * Ks.KP1 * Ks.KP2 * Ks.KP3 - H^3
    PhosDenom = H^3 + Ks.KP1 * H^2 + Ks.KP1 * Ks.KP2 * H + Ks.KP1 * Ks.KP2 * Ks.KP3
    PAlk = PT * PhosNum / PhosDenom
    SiAlk = SiT * Ks.KSi / (Ks.KSi + H)
    Alk_H2S = H2ST * Ks.KH2S / (Ks.KH2S + H)
    Alk_NH3 = NH4T * Ks.KNH3 / (Ks.KNH3 + H)
    Hfree = H / (1 + ST / Ks.KS)
    HSO₄ = ST / (1 + Ks.KS / Hfree)
    HF = FT / (1 + Ks.KF / Hfree)

    return TA - CAlk - BAlk - OH - PAlk - SiAlk - Alk_H2S - Alk_NH3 + Hfree + HSO₄ + HF
end


function pH_from_TA_DIC(TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # 1. Create a temporary function where pH is the only input
    f(pH) = solve_pH_from_TA_DIC(pH, TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    
    # 2. Use ForwardDiff to automatically generate the exact derivative function
    df(pH) = ForwardDiff.derivative(f, pH)
    
    # 3. Pass both the function and its derivative to find_zero, using Newton's method
    initial_guess = 8.0 + zero(TA) + zero(DIC)
    return find_zero((f, df), initial_guess, Roots.Newton())
end


## Equation 1.1.9
"""
Calculating CO₂ from H⁺ and DIC
Equation 1.1.9 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_CO₂(H, DIC, Ks)
    return DIC / (1 + Ks.K1 / H + Ks.K1 * Ks.K2 / H^2)
end


## Equation 1.1.10
"""
Calculating HCO₃ from H⁺ and DIC
Equation 1.1.10 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_HCO₃(H, DIC, Ks)
    return DIC / (1 + H / Ks.K1 + Ks.K2 / H)
end


## Equation 1.1.11
"""
Calculating CO₃ from H⁺ and DIC
Equation 1.1.11 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_CO₃(H, DIC, Ks)
    return DIC / (1 + H / Ks.K2 + H^2 / (Ks.K1 * Ks.K2))
end


## Equation 1.5.80
"""
Calculating TA
Equation 1.5.80 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_TA(H, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks; mode="multi")
    Denom = H^2 + Ks.K1 * H + Ks.K1 * Ks.K2
    CAlk = DIC * Ks.K1 * (H + 2 * Ks.K2) / Denom
    BAlk = BT * Ks.KB / (Ks.KB + H)
    OH = Ks.KW / H
    PhosNum = Ks.KP1 * Ks.KP2 * H + 2 * Ks.KP1 * Ks.KP2 * Ks.KP3 - H^3
    PhosDenom = H^3 + Ks.KP1 * H^2 + Ks.KP1 * Ks.KP2 * H + Ks.KP1 * Ks.KP2 * Ks.KP3
    PAlk = PT * PhosNum / PhosDenom
    SiAlk = SiT * Ks.KSi / (Ks.KSi + H)
    Alk_H2S = H2ST * Ks.KH2S / (Ks.KH2S + H)
    Alk_NH3 = NH4T * Ks.KNH3 / (Ks.KNH3 + H)
    Hfree = H / (1 + ST / Ks.KS)
    HSO₄ = ST / (1 + Ks.KS / Hfree)
    HF = FT / (1 + Ks.KF / Hfree)

    TA = CAlk + BAlk + OH + PAlk + SiAlk + Alk_H2S + Alk_NH3 - Hfree - HSO₄ - HF

    if mode == "multi"
        return TA, CAlk, BAlk, PAlk, SiAlk, OH, Hfree, HSO₄, HF, Alk_H2S, Alk_NH3
    else
        return TA
    end
end


# Equation C.4.14
"""
Calculating CO₂ from fugacity
Equation C.4.14 from Zeebe & Wolf-Gladrow, 2001, Appendix C
"""
function fCO₂_to_CO₂(fCO₂, Ks)
    return fCO₂ * Ks.K0
end 


# Equation C.4.14
"""
Calculating fugacity from CO₂
Equation C.4.14 from Zeebe & Wolf-Gladrow, 2001, Appendix C
"""
function CO₂_to_fCO₂(CO₂, Ks)
    return CO₂ / Ks.K0
end 


# pCO₂ --> fCO₂
"""
Calculating fCO₂ from pCO₂
Taken from MatLab CO2SYS

Assumes a pressure of or near 1 atm, otherwise, the exponential pressure term
will impact calculations (Weiss, R. F., Marine Chemistry 2:203-215, 1974).

Intended for a mixture of CO₂ and air at 1 atm (low CO₂ concentrations).

Δ & B are in cm³/mol
"""
function pCO₂_to_fCO₂(pCO₂, T)
    Tₖ = T + 273.15
    P = 1.01325 # in bar
    RT = 83.14472 * Tₖ # originally used R = 83.14462618, however switched to 83.14472 to match CO2SYS
    a₀, a₁, a₂, a₃ = (-1636.75, 12.0408, -3.27957e-2, 3.16528e-05)
    b₀, b₁ = (57.7, -0.118)
    B = a₀ + a₁ * Tₖ+ a₂ * Tₖ^2 + a₃ * Tₖ^3
    Δ = b₀ + b₁ * Tₖ
    return pCO₂ * exp(P * (B + 2 * Δ)/RT)
end


# fCO₂ --> pCO₂
"""
Calculating pCO₂ from fCO₂
Taken from MatLab CO2SYS

Assumes a pressure of or near 1 atm, otherwise, the exponential pressure term
will impact calculations (Weiss, R. F., Marine Chemistry 2:203-215, 1974).

Intended for a mixture of CO₂ and air at 1 atm (low CO₂ concentrations).

Δ & B are in cm³/mol
"""
function fCO₂_to_pCO₂(fCO₂, T)
    Tₖ = T + 273.15
    P = 1.01325 # in bar
    RT = 83.14472 * Tₖ # originally used R = 83.14462618, however switched to 83.14472 to match CO2SYS
    a₀, a₁, a₂, a₃ = (-1636.75, 12.0408, -3.27957e-2, 3.16528e-05)
    b₀, b₁ = (57.7, -0.118)
    B = a₀ + a₁ * Tₖ + a₂ * Tₖ^2 + a₃ * Tₖ^3
    Δ = b₀ + b₁ * Tₖ
    return fCO₂ / exp(P * (B + 2 * Δ) / RT)
end


# The Carbon Calculator
"""
Calculates carbon system from any two of the following: 
CO₂, HCO₃⁻, CO₃²⁻, DIC, TA, pH
Returns everything on the total scale.
"""
function C_calculator(;
    pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing, 
    CO₃=nothing, fCO₂=nothing, pCO₂=nothing, fH=nothing, BT=0.0, PT=0.0, SiT=0.0,
    ST=0.0, FT=0.0, H2ST=0.0, NH4T=0.0, Ks=nothing, T_in=25.0, S_in=35.0, 
    kwargs...)

    # If fCO₂ is given but CO₂ is not, calculate CO₂:
    if isnothing(CO₂)
        if !isnothing(fCO₂)
            CO₂ = fCO₂_to_CO₂(fCO₂, Ks)
        elseif !isnothing(pCO₂)
            # Calculate fCO2 once and store it!
            fCO₂ = pCO₂_to_fCO₂(pCO₂, T_in)
            CO₂ = fCO₂_to_CO₂(fCO₂, Ks)
        end
    end

    # Calculations based on logic in Zeebe & Wolf-Gladrow, 2001, Appendix B

    # 1. CO₂ and pH given; calculate H⁺ and DIC
    if !isnothing(CO₂) && !isnothing(pHtot)
        H = 10.0^(-pHtot)
        DIC = DIC_from_CO₂_pH(CO₂, pHtot, Ks)
    # 2. CO₂ and HCO₃ given; calculate H⁺, pHtot, and DIC
    elseif !isnothing(CO₂) && !isnothing(HCO₃)
        H = H_from_CO₂_HCO₃(CO₂, HCO₃, Ks)
        pHtot = -log10(H)
        DIC = DIC_from_CO₂_pH(CO₂, pHtot, Ks)
    # 3. CO₂ and CO₃; calculate H⁺ and DIC
    elseif !isnothing(CO₂) && !isnothing(CO₃)
        H = H_from_CO₂_CO₃(CO₂, CO₃, Ks)
        DIC = DIC_from_CO₂_pH(CO₂, -log10(H), Ks)
    # 4. CO₂ and TA; calculate H⁺, pHtot, and DIC
    elseif !isnothing(CO₂) && !isnothing(TA)
        pHtot = pH_from_CO₂_TA(CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
        H = 10.0^(-pHtot)
        DIC = DIC_from_CO₂_pH(CO₂, pHtot,Ks)
    # 5. CO₂ and DIC; calculate H⁺
    elseif !isnothing(CO₂) && !isnothing(DIC)
        H = H_from_CO₂_DIC(CO₂, DIC, Ks)
    #6. pH and HCO₃; calculate H⁺ and DIC
    elseif !isnothing(pHtot) && !isnothing(HCO₃)
        H = 10.0^(-pHtot)
        DIC = DIC_from_pH_HCO₃(pHtot, HCO₃, Ks)
    #7. pH and CO₃; calculate H⁺ and DIC
    elseif !isnothing(pHtot) && !isnothing(CO₃)
        H = 10.0^(-pHtot)
        DIC = DIC_from_pH_CO₃(pHtot, CO₃, Ks)
    # 8. pH and TA; calculate H⁺ and DIC
    elseif !isnothing(pHtot) && !isnothing(TA)
        H = 10.0^(-pHtot)
        DIC = DIC_from_pH_TA(pHtot, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # 9. pH and DIC; calculate H⁺
    elseif !isnothing(pHtot) && !isnothing(DIC)
        H = 10.0^(-pHtot)
    # 10. HCO₃ and CO₃; calculate H⁺ and DIC
    elseif !isnothing(HCO₃) && !isnothing(CO₃)
        H = H_from_HCO₃_CO₃(HCO₃, CO₃, Ks)
        DIC = DIC_from_pH_CO₃(-log10(H), CO₃, Ks)
    # 11. HCO₃ and TA; calculate H⁺ and DIC
    elseif !isnothing(HCO₃) && !isnothing(TA)
        H = H_from_HCO₃_TA(HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
        DIC = HCO₃ * (H / Ks.K1 + 1.0 + Ks.K2 / H) 
        pHtot = -log10(H)
    # 12. HCO₃ and DIC; calculate H⁺ and pHtot
    elseif !isnothing(HCO₃) && !isnothing(DIC)
        pHtot = pH_from_HCO₃_DIC(HCO₃,DIC, Ks)
        H = 10.0^(-pHtot)
    #13. CO₃ and TA; calculate H⁺ and DIC
    elseif !isnothing(CO₃) && !isnothing(TA)
        H = H_from_CO₃_TA(CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
        DIC = CO₃ * (H^2 / (Ks.K1 * Ks.K2) + H / Ks.K2 + 1.0)
        pHtot = -log10(H)
    #14. CO₃ and DIC; calculate H⁺
    elseif !isnothing(CO₃) && !isnothing(DIC)
        H = H_from_CO₃_DIC(CO₃, DIC, Ks)
    #15. TA and DIC; calculate H⁺ and pHtot
    elseif !isnothing(TA) && !isnothing(DIC)
        pHtot = pH_from_TA_DIC(TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
        H = 10.0^(-pHtot)
    end

    if isnothing(CO₂)
        CO₂ = calc_CO₂(H, DIC, Ks)
    end
    if isnothing(fCO₂)
        fCO₂ = CO₂_to_fCO₂(CO₂, Ks)
    end
    if isnothing(pCO₂)
        pCO₂ = fCO₂_to_pCO₂(fCO₂, T_in)
    end
    if isnothing(HCO₃)
        HCO₃ = calc_HCO₃(H, DIC, Ks)
    end
    if isnothing(CO₃)
        CO₃ = calc_CO₃(H, DIC, Ks)
    end
    
    (TA, CAlk, BAlk, PAlk, SiAlk, OH, Hfree, HSO₄, HF, Alk_H2S, Alk_NH3) = calc_TA(
        H, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks; mode="multi"
    )

    if isnothing(pHtot)
        pHtot = -log10(H)
    end

    # Return only core variables on the Total scale
    return (; 
        pHtot=pHtot, TA=TA, DIC=DIC, CO₂=CO₂, H=H, HCO₃=HCO₃, CO₃=CO₃,
        fCO₂=fCO₂, pCO₂=pCO₂, CAlk=CAlk, BAlk=BAlk, PAlk=PAlk, SiAlk=SiAlk, OH=OH,
        Hfree=Hfree, HSO₄=HSO₄, HF=HF, Alk_H2S=Alk_H2S, Alk_NH3=Alk_NH3
    )
end

# Calculate the Revelle Factor
"""
Calculating the Revelle Factor from CO₂ and DIC:
ΔpCO₂ / ΔDIC
"""
function calc_revelle_factor(TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # 1. First, determine baseline fCO2
    pH_base = pH_from_TA_DIC(TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    fCO₂_base = calc_CO₂(10.0^(-pH_base), DIC, Ks) / Ks.K0
    
    # 2. Create a function that calculates fCO2 entirely from a given DIC
    function fCO₂_from_DIC(d)
        pH_temp = pH_from_TA_DIC(TA, d, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
        return calc_CO₂(10.0^(-pH_temp), d, Ks) / Ks.K0
    end
    
    # 3. Get the EXACT rate of change of fCO2 with respect to DIC
    dfCO₂_dDIC = ForwardDiff.derivative(fCO₂_from_DIC, DIC)
    
    # 4. Return the Revelle fraction (∂fCO2 / ∂DIC) * (DIC / fCO2)
    return dfCO₂_dDIC * (DIC / fCO₂_base)
end

"""
Calculates the TA Buffer Capacity (∂TA / ∂pH) using Automatic Differentiation
"""
function calc_buffer_capacity(pH, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # 1. Create a temporary function where pH is the ONLY input
    # Note: Pass "single" so it only returns the TA value, not the breakdown
    f_TA(p) = calc_TA(10.0^(-p), DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks; mode="multi")

    
    # 2. Get the exact derivative!
    dTA_dpH = ForwardDiff.derivative(f_TA, pH)
    
    return dTA_dpH
end

end # module

