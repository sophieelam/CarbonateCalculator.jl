module Carbon
using Roots
using ForwardDiff

export C_calculator, fCO₂_to_CO₂,
CO₂_to_fCO₂, fCO₂_to_pCO₂, pCO₂_to_fCO₂

# `calc_revelle_factor` and `calc_buffer_capacity` are deliberately absent from the export list. Both
# are superseded by `calc_gradient`, which reproduces them exactly, and both are easy to misuse
# from outside: they run before `_rescale_to_unit`, so they take mol/kg, where the obvious thing
# for a caller to do is pass values off `result.val` — which are in the reporting unit, and give
# a plausible answer about 1% wrong. Kept in place, and still reachable within the package as
# `Carbon.calc_revelle_factor`.

# One knob for the iterative solves below. `const` matters: an untyped global would make
# every `ROOT_METHOD()` a dynamic dispatch on the hot path.
const ROOT_METHOD = Roots.Newton

"Where the iterative solves start. Seawater sits near here, so most inputs converge at once."
const DEFAULT_PH_GUESS = 8.0

"The pH range searched when `DEFAULT_PH_GUESS` fails. Wider than any solvable input."
const PH_SEARCH_RANGE = (2.0, 13.0)

"Resolution of that search, in pH units. Fine enough to leave Newton well inside the basin."
const PH_SEARCH_STEP = 0.5

"""
    _newton_state_type(Ks, concentrations...)

The type the Newton solver's internal state has to hold.

Promotes over the equilibrium constants as well as the concentrations, because either can
carry `ForwardDiff.Dual`s. An uncertainty on `temp_c` or `sal` reaches the answer only through
the constants, leaving the concentrations `Float64`; promoting over the concentrations alone
would then start the solve at `Float64` against a `Dual` residual, and `Roots` raises
`MethodError: no method matching Float64(::Dual)` from inside `init_state`.

Resolved at compile time: every argument's type is known, so this folds to a constant.
"""
_newton_state_type(Ks::NamedTuple, concentrations...) =
    promote_type(map(typeof, concentrations)..., map(typeof, Tuple(values(Ks)))...)

"""
    _bracketed_guess(residual, T) -> pH

A starting point inside the root's basin, found by scanning [`PH_SEARCH_RANGE`](@ref) for a
sign change and taking the midpoint of the interval that has one.

Reached only when the default start fails, so the scan costs nothing on inputs that converge
without it.

Comparisons work directly on `ForwardDiff.Dual`s — Julia orders them by value — so nothing has
to be stripped. The result is converted to `T` so a `Float64` start never meets a `Dual`
residual (see [`_newton_state_type`](@ref)); its partials are irrelevant, since Newton
recovers those from the residual itself.

Falls back to the default guess when no sign change is found, leaving the caller to fail as
it otherwise would.
"""
function _bracketed_guess(residual, ::Type{T}) where {T}
    low, high = PH_SEARCH_RANGE
    steps = round(Int, (high - low) / PH_SEARCH_STEP)

    left = convert(T, low)
    f_left = residual(left)

    for step in 1:steps
        right = convert(T, low + step * PH_SEARCH_STEP)
        f_right = residual(right)
        (f_left > 0) == (f_right > 0) || return (left + right) / 2
        left, f_left = right, f_right
    end

    return convert(T, DEFAULT_PH_GUESS)
end

"""
    _solve_pH(residual, derivative, T) -> pH

Newton from [`DEFAULT_PH_GUESS`](@ref), falling back to a bracketed start when that diverges.

pH 8 converges immediately for seawater but diverges when the answer is far from it — beyond
roughly pH 10.5 for `TA+DIC` and 9.5 for `CO₂+TA`. The first attempt goes through
`Roots.solve`, which reports failure as `NaN` rather than throwing, so the retry is ordinary
control flow; an input with no root still raises from the second attempt.

**Newton does the converging in both cases, deliberately.** `find_zero` given an *interval*
returns a plain `Float64` whatever it is handed, so a `Dual` passed in comes back stripped of
its partials and any uncertainty propagated through this path arrives as exactly zero.
Newton's iteration is ordinary arithmetic, so `Dual`s survive it — bracketing places the
starting point, and never finds the root.

Only valid for residuals that are monotonic in pH, which is every pair that reaches here.
"""
function _solve_pH(residual, derivative, ::Type{T}) where {T}
    from_default = Roots.solve(
        Roots.ZeroProblem((residual, derivative), convert(T, DEFAULT_PH_GUESS)),
        ROOT_METHOD()
    )
    isnan(from_default) || return from_default

    return find_zero((residual, derivative), _bracketed_guess(residual, T), ROOT_METHOD())
end

"""
    _quadratic_roots(a, b, c) -> (root, root)

Both roots of `a*x^2 + b*x + c`, in no particular order.

Uses the cancellation-avoiding form rather than the textbook quadratic formula. One of
`-b ± sqrt(b^2 - 4ac)` always subtracts two nearly equal numbers when `b^2 ≫ 4ac`, losing
precision in exactly one of the two roots. Forming `q` with the signs aligned and taking
`q/a` and `c/q` makes both of them quotients of well-conditioned quantities.

The carbonate pairs that reduce to a quadratic use this instead of a root-finder, which keeps them faster and differentiable.
Callers pick the root they want by sign or magnitude; which of the pair is which depends on
`sign(b)`, so do not rely on the order.
"""
function _quadratic_roots(a, b, c)
    q = -(b + sign(b) * sqrt(b^2 - 4a * c)) / 2
    return (q / a, c / q)
end

"""
    _positive_root(a, b, c)

The positive root of `a*H^2 + b*H + c`, for the case where `a > 0` and `c < 0`.

Those signs put the product of the roots (`c/a`) below zero, so the two straddle zero and
exactly one is physical — which is why `max` is enough to identify it. Callers whose
constant term is *positive* have two positive roots and no such shortcut; they take
`_quadratic_roots` directly and choose on magnitude.
"""
_positive_root(a, b, c) = max(_quadratic_roots(a, b, c)...)


"""
#1: Calculating DIC from CO₂ and pH
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""

function DIC_from_CO₂_pH(CO₂, pH, Ks)
    H = 10.0^(-pH)
    return CO₂ * (1.0 + Ks.K1/H + Ks.K1*Ks.K2/H^2)
end


"""
#2: Calculating H⁺ from CO₂ and HCO₃⁻
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function H_from_CO₂_HCO₃(CO₂, HCO₃, Ks)
    return (Ks.K1 * CO₂) / HCO₃
end


"""
#3: Calculating H⁺ from CO₂ and CO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function H_from_CO₂_CO₃(CO₂, CO₃, Ks)
    # Using abs() just as a safety net against tiny floating-point noise around zero
    return sqrt(abs((Ks.K1 * Ks.K2 * CO₂) / CO₃))
end


"""
    _non_carbonate_alkalinity(H, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

Every contribution to alkalinity except the carbonate one, at a given [H⁺]: borate, hydroxide,
phosphate, silicate, sulphide, ammonia, and the free-proton terms that subtract.

Shared by the four places alkalinity is assembled — the residuals for `CO₂+TA` and `TA+DIC`,
the `DIC`-from-`pH`-and-`TA` rearrangement, and `calc_TA_components` — which agree on all of
these and differ only in where CAlk comes from and which way the sum runs. Changing how
borate or phosphate is handled therefore means changing it here alone.

Returned as a NamedTuple because `calc_TA_components` reports the breakdown, not just the sum.
"""
function _non_carbonate_alkalinity(H, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    PhosNum = Ks.KP1 * Ks.KP2 * H + 2 * Ks.KP1 * Ks.KP2 * Ks.KP3 - H^3
    PhosDenom = H^3 + Ks.KP1 * H^2 + Ks.KP1 * Ks.KP2 * H + Ks.KP1 * Ks.KP2 * Ks.KP3
    Hfree = H / (1 + ST / Ks.KS)

    return (BAlk = BT * Ks.KB / (Ks.KB + H),
            OH = Ks.KW / H,
            PAlk = PT * PhosNum / PhosDenom,
            SiAlk = SiT * Ks.KSi / (Ks.KSi + H),
            Alk_H2S = H2ST * Ks.KH2S / (Ks.KH2S + H),
            Alk_NH3 = NH4T * Ks.KNH3 / (Ks.KNH3 + H),
            Hfree = Hfree,
            HSO₄ = ST / (1 + Ks.KS / Hfree),
            HF = FT / (1 + Ks.KF / Hfree))
end

"""
    _alkalinity_residual(TA, CAlk, contributions)

What is left of `TA` once carbonate and every other contribution is taken out. Zero at the
solution, which is what the iterative solves search for.

Passing `CAlk = 0` gives the carbonate alkalinity implied by a `TA` — the
rearrangement `DIC_from_pH_TA` needs.
"""
_alkalinity_residual(TA, CAlk, a) =
    TA - CAlk - a.BAlk - a.OH - a.PAlk - a.SiAlk - a.Alk_H2S - a.Alk_NH3 +
    a.Hfree + a.HSO₄ + a.HF

"Carbonate alkalinity, HCO₃⁻ + 2CO₃²⁻, from DIC and [H⁺]."
_carbonate_alkalinity(H, DIC, Ks) =
    DIC * Ks.K1 * (H + 2 * Ks.K2) / (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2)


"""
#4: pH from CO₂ and TA

Alkalinity as CO₂ and this pH imply it, less the alkalinity actually supplied. Solved by
Newton with an analytic `ForwardDiff` derivative; see [`_solve_pH`](@ref).
"""
function solve_pH_from_CO₂_TA(pH, CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

    H = 10.0^(-pH)

    # CAlk comes from CO₂ here rather than from DIC, which is the only thing distinguishing
    # this residual from `solve_pH_from_TA_DIC`. The round trip through K0 is CO2SYS's and is
    # kept because dropping it moves the last bits.
    fCO₂ = CO₂ / Ks.K0
    HCO₃ = Ks.K0 * Ks.K1 * fCO₂ / H
    CO₃ = Ks.K0 * Ks.K1 * Ks.K2 * fCO₂ / H^2
    CAlk = HCO₃ + 2 * CO₃

    return _alkalinity_residual(TA, CAlk,
        _non_carbonate_alkalinity(H, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks))
end

function pH_from_CO₂_TA(CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    T = _newton_state_type(Ks, CO₂, TA, BT)

    f(pH) = solve_pH_from_CO₂_TA(pH, CO₂, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    df(pH) = ForwardDiff.derivative(f, pH)

    return _solve_pH(f, df, T)
end


"""
#5: Calculating H⁺ from CO₂ and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B

`DIC·H² = CO₂·(H² + K₁H + K₁K₂)` rearranges to `(DIC − CO₂)H² − CO₂K₁H − CO₂K₁K₂ = 0`, so
the answer is a quadratic root. `DIC > CO₂` always, which makes the
leading coefficient positive and the constant negative, so one positive root.

Solved in closed form because iterating on this residual does not work: it is of order 1e-19
near the root, so any absolute convergence test is satisfied before the iteration has done
anything.
"""
function H_from_CO₂_DIC(CO₂, DIC, Ks) 
    return _positive_root(
        DIC - CO₂, 
        -CO₂ * Ks.K1, 
        -CO₂ * Ks.K1 * Ks.K2
    )
end


"""
#6: Calculating DIC from pH and HCO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function DIC_from_pH_HCO₃(pH, HCO₃, Ks)
    H = 10.0^(-pH)
    return HCO₃ * (1.0 + H / Ks.K1 + Ks.K2 / H)
end


"""
#7: Calculating DIC from pH and CO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function DIC_from_pH_CO₃(pH, CO₃, Ks)
    H = 10.0^(-pH)
    return CO₃ * (1 + H / Ks.K2 + H^2 / (Ks.K1 * Ks.K2))
end


"""
#8: Calculating DIC from pH and TA
Taken from MatLab CO2SYS
"""
function DIC_from_pH_TA(pH, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    H = 10.0^(-pH)

    # Carbonate alkalinity is whatever is left of TA once everything else is accounted for,
    # which is the residual with no carbonate term subtracted.
    CAlk = _alkalinity_residual(TA, zero(TA),
        _non_carbonate_alkalinity(H, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks))

    # `_carbonate_alkalinity` inverted for DIC.
    return CAlk * (H^2 + Ks.K1 * H + Ks.K1 * Ks.K2) / (Ks.K1 * (H + 2.0 * Ks.K2))
end


"""
#9: Calculating CO₂ from pH and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function CO₂_from_pH_DIC(pH, DIC, Ks)
    H = 10.0^(-pH)
    return DIC / (1.0 + Ks.K1 / H + Ks.K1 * Ks.K2 / H^2)
end


"""
#10: Calculating H⁺ from HCO₃ and CO₃
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function H_from_HCO₃_CO₃(HCO₃, CO₃, Ks) 
    return Ks.K2 * HCO₃ / CO₃
end


"""
#11: Calculating H⁺ from HCO₃ and TA
Zeebe & Wolf-Gladrow, 2001, Appendix B
"""
function solve_H_from_HCO₃_TA(H, HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    temp_DIC = HCO₃ * (H / Ks.K1 + 1.0 + Ks.K2 / H)
    calc_TA_val = calc_TA(H, temp_DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    return calc_TA_val - TA
end

function H_from_HCO₃_TA(HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    T = _newton_state_type(Ks, HCO₃, TA, BT)

    f(pH) = solve_H_from_HCO₃_TA(10.0^(-pH), HCO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    df(pH) = ForwardDiff.derivative(f, pH)

    return 10.0^(-_solve_pH(f, df, T))
end

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

    # Unlike the other quadratics here the constant term is *positive*, so both roots are
    # positive and `_positive_root` does not apply — there is no sign to choose on. HCO₃
    # close to DIC leaves no real root at all, which is a physically inconsistent pair
    # rather than a failure of the algebra.
    b^2 - 4a * c < 0 && return NaN

    # The smaller root, i.e. the higher pH. Taken from `_quadratic_roots` rather than
    # written out because the smaller one is precisely the root the textbook formula
    # computes badly: it comes from -b - sqrt(disc), a subtraction of two close numbers
    # once 4ac is small next to b². At HCO₃ = 400, DIC = 2400 µmol that costs four digits.
    smaller_H = min(_quadratic_roots(a, b, c)...)

    return -log10(smaller_H) # This doesn't match CBsyst, but matches test suite
end


"""
#13: Calculating H⁺ from CO₃ and TA
Zeebe & Wolf-Gladrow, 2001, Appendix B
Uses Roots.Brent() (same as CBsys) to circumnavigate bracketing issues with
root finding. However, this only works for pH values 5 < pH < 10.
"""
function solve_H_from_CO₃_TA(H, CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    temp_DIC = CO₃ * (H^2 / (Ks.K1 * Ks.K2) + H / Ks.K2 + 1)
    calc_TA_val = calc_TA(H, temp_DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    return calc_TA_val - TA
end

function H_from_CO₃_TA(CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    T = _newton_state_type(Ks, CO₃, TA, BT)

    f(pH) = solve_H_from_CO₃_TA(10.0^(-pH), CO₃, TA, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    df(pH) = ForwardDiff.derivative(f, pH)

    # Not `_solve_pH`, which assumes a single root. This residual is not monotonic in pH: at
    # fixed CO₃ the implied DIC grows as H² while the free-proton terms eventually pull TA
    # back down, so every input has two solutions and both are physically valid. Which one
    # Newton reaches depends on where it starts, and this pair carries nothing to choose on —
    # answers should not be trusted.

    # TODO: add a warning about this when this pair is used.
    initial_guess = convert(T, DEFAULT_PH_GUESS)
    sol_pH = find_zero((f, df), initial_guess, ROOT_METHOD())

    return 10.0^(-sol_pH)
end

"""
#14: Calculating H⁺ from CO₃ and DIC
Zeebe & Wolf-Gladrow, 2001, Appendix B

`CO₃·(1 + H/K₂ + H²/(K₁K₂)) = DIC` is a quadratic in H with coefficients
`CO₃/(K₁K₂)`, `CO₃/K₂` and `CO₃ − DIC`. `DIC > CO₃` always, so the constant term is negative
and the positive root is the physical one.
"""
function H_from_CO₃_DIC(CO₃, DIC, Ks)
    return _positive_root(
        CO₃ / (Ks.K1 * Ks.K2), 
        CO₃ / Ks.K2, 
        CO₃ - DIC
    )
end


"""
#15: pH from TA and DIC

The most-used pair. Alkalinity as DIC and this pH imply it, less the alkalinity supplied.
Monotonic in pH, so it has exactly one root.
"""
function solve_pH_from_TA_DIC(pH, TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    H = 10.0^(-pH)

    return _alkalinity_residual(TA, _carbonate_alkalinity(H, DIC, Ks),
        _non_carbonate_alkalinity(H, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks))
end


function pH_from_TA_DIC(TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    T = _newton_state_type(Ks, TA, DIC, BT)

    f(pH) = solve_pH_from_TA_DIC(pH, TA, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    df(pH) = ForwardDiff.derivative(f, pH)

    return _solve_pH(f, df, T)
end


"""
Calculating CO₂ from H⁺ and DIC
Equation 1.1.9 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_CO₂(H, DIC, Ks)
    return DIC / (1 + Ks.K1 / H + Ks.K1 * Ks.K2 / H^2)
end


"""
Calculating HCO₃ from H⁺ and DIC
Equation 1.1.10 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_HCO₃(H, DIC, Ks)
    return DIC / (1 + H / Ks.K1 + Ks.K2 / H)
end


"""
Calculating CO₃ from H⁺ and DIC
Equation 1.1.11 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_CO₃(H, DIC, Ks)
    return DIC / (1 + H / Ks.K2 + H^2 / (Ks.K1 * Ks.K2))
end


"""
Calculating TA components
Equation 1.5.80 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_TA_components(H, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    CAlk = _carbonate_alkalinity(H, DIC, Ks)
    (; BAlk, OH, PAlk, SiAlk, Alk_H2S, Alk_NH3, Hfree, HSO₄, HF) =
        _non_carbonate_alkalinity(H, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

    TA = CAlk + BAlk + OH + PAlk + SiAlk + Alk_H2S + Alk_NH3 - Hfree - HSO₄ - HF

    return TA, CAlk, BAlk, PAlk, SiAlk, OH, Hfree, HSO₄, HF, Alk_H2S, Alk_NH3
end

"""
Calculating TA
Equation 1.5.80 from Zeebe & Wolf-Gladrow, 2001, Chapter 1
"""
function calc_TA(H, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    (TA, _, _, _, _, _, _, _, _, _, _) = calc_TA_components(H, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

    return TA
end


"""
Calculating CO₂ from fugacity
Equation C.4.14 from Zeebe & Wolf-Gladrow, 2001, Appendix C
"""
function fCO₂_to_CO₂(fCO₂, Ks)
    return fCO₂ * Ks.K0
end 


"""
Calculating fugacity from CO₂
Equation C.4.14 from Zeebe & Wolf-Gladrow, 2001, Appendix C
"""
function CO₂_to_fCO₂(CO₂, Ks)
    return CO₂ / Ks.K0
end 


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
    RT = 83.14472 * Tₖ # 83.14472, not the DOEv2 83.14462618, to match CO2SYS
    a₀, a₁, a₂, a₃ = (-1636.75, 12.0408, -3.27957e-2, 3.16528e-05)
    b₀, b₁ = (57.7, -0.118)
    B = a₀ + a₁ * Tₖ+ a₂ * Tₖ^2 + a₃ * Tₖ^3
    Δ = b₀ + b₁ * Tₖ
    return pCO₂ * exp(P * (B + 2 * Δ)/RT)
end


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
    RT = 83.14472 * Tₖ # 83.14472, not the DOEv2 83.14462618, to match CO2SYS
    a₀, a₁, a₂, a₃ = (-1636.75, 12.0408, -3.27957e-2, 3.16528e-05)
    b₀, b₁ = (57.7, -0.118)
    B = a₀ + a₁ * Tₖ + a₂ * Tₖ^2 + a₃ * Tₖ^3
    Δ = b₀ + b₁ * Tₖ
    return fCO₂ / exp(P * (B + 2 * Δ) / RT)
end


"""
What to tell someone whose input reaches the carbonate solver without determining it.

`whole_system` requires two constraints overall rather than two *carbonate* ones, so
`whole_system(pHtot = 8.1, δBOH₄ = 16.0)` determines boron and its isotopes while leaving
carbon with pH alone. No branch of the solver matches, and the message has to name the
missing input rather than let an internal variable go undefined.
"""
function _underdetermined_carbon(pHtot, DIC, TA, CO₂, HCO₃, CO₃)
    supplied = [name for (name, value) in (("pH", pHtot), ("DIC", DIC), ("TA", TA),
                                           ("CO₂", CO₂), ("HCO₃", HCO₃), ("CO₃", CO₃))
                if !isnothing(value)]

    return "the carbonate system needs two of TA, DIC, pH, CO₂ (or pCO₂/fCO₂), HCO₃, CO₃, " *
           "but was given $(isempty(supplied) ? "none" : join(supplied, " and ")).\n" *
           "A scope that includes :carbon computes the carbonate system, so it has to be " *
           "determined too. To solve boron and its isotopes alone, use boron_system or " *
           "boron_isotopes, whose scope leaves carbon out."
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
    ST=0.0, FT=0.0, H2ST=0.0, NH4T=0.0, Ks=nothing, temp_c=25.0, sal=35.0, 
    kwargs...)

    # If fCO₂ is given but CO₂ is not, calculate CO₂:
    if isnothing(CO₂)
        if !isnothing(fCO₂)
            CO₂ = fCO₂_to_CO₂(fCO₂, Ks)
        elseif !isnothing(pCO₂)
            # Calculate fCO2 once and store it!
            fCO₂ = pCO₂_to_fCO₂(pCO₂, temp_c)
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
    else
        throw(ArgumentError(_underdetermined_carbon(pHtot, DIC, TA, CO₂, HCO₃, CO₃)))
    end

    if isnothing(CO₂)
        CO₂ = calc_CO₂(H, DIC, Ks)
    end
    if isnothing(fCO₂)
        fCO₂ = CO₂_to_fCO₂(CO₂, Ks)
    end
    if isnothing(pCO₂)
        pCO₂ = fCO₂_to_pCO₂(fCO₂, temp_c)
    end
    if isnothing(HCO₃)
        HCO₃ = calc_HCO₃(H, DIC, Ks)
    end
    if isnothing(CO₃)
        CO₃ = calc_CO₃(H, DIC, Ks)
    end
    
    (TA, CAlk, BAlk, PAlk, SiAlk, OH, Hfree, HSO₄, HF, Alk_H2S, Alk_NH3) = calc_TA_components(
        H, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks
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
    calc_buffer_capacity(pH, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

Buffer capacity ∂TA/∂pH, in mol/kg per pH unit, by automatic differentiation.
"""
function calc_buffer_capacity(pH, DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)
    # `calc_TA`, not `calc_TA_components`: ForwardDiff needs a scalar-valued function, and the
    # eleven-value breakdown gives `MethodError: no method matching extract_derivative`.
    f_TA(p) = calc_TA(10.0^(-p), DIC, BT, PT, SiT, ST, FT, H2ST, NH4T, Ks)

    return ForwardDiff.derivative(f_TA, pH)
end

end # module

