# How a system is solved: the stages of the calculation, in the order they run.
#
# The tables these read — which parameters exist, which constrain which subsystem, which are
# concentrations — live in parameters.jl. This file is the procedure; that file describes the
# inputs.
#
# Each stage is a named function rather than a section of one long block, so the order of the
# calculation is visible in `_solve_core` and each step can be read on its own.

"Groups in `inputs` that are supplied, for one subsystem."
function _constraints(inputs, subsystem::Symbol)
    n = 0
    for (group, members) in pairs(PARAMETER_GROUPS)
        getproperty(GROUP_SUBSYSTEM, group) === subsystem || continue
        any(m -> haskey(inputs, m) && !isnothing(inputs[m]), members) && (n += 1)
    end
    return n
end

"""
    _first_solvable(ps) -> :boron | :isotopes | :carbon | :none

Which subsystem to solve first.

The three calculators are always all run — the choice is only what order, because pH is the
shared unknown. Whichever subsystem has enough constraints to pin pH goes first, and the
other two then follow with pH already known.

Carbon needs two constraints; boron and isotopes need only one, because their totals (`BT`,
`δBT`) always have a value, so a single speciated measurement closes them. A known pH
satisfies everything at once.
"""
function _first_solvable(ps)
    # If pH is *known*, return early: all three subsystems are then immediately solvable,
    # the order is irrelevant, and boron leads only because something has to.
    isnothing(ps.pHtot) || return :boron

    # Otherwise take whichever subsystem the inputs determine. At most one of these can
    # hold: more than two constraints is over-determined, and the wrappers reject that
    # before the core is reached, so the order between them is not a precedence rule.
    _constraints(ps, :boron) >= CONSTRAINTS_NEEDED.boron && return :boron
    _constraints(ps, :isotopes) >= CONSTRAINTS_NEEDED.isotopes && return :isotopes
    # The common case reaches here: two carbonate parameters and no boron speciation, so
    # carbon solves first and yields pH, and boron and the isotopes follow from it.
    _constraints(ps, :carbon) >= CONSTRAINTS_NEEDED.carbon && return :carbon

    return :none
end

"What to tell someone whose input does not determine the system, built from the table."
function _underdetermined_message()
    routes = [
        "pH on any scale",
        "two of $(join(("TA", "DIC", "CO₂ (or pCO₂/fCO₂)", "HCO₃", "CO₃ (or ΩA/ΩC)"), ", "))",
        "BT with one of BOH₃, BOH₄",
        "δBT with one of δBOH₃, δBOH₄ (or the A equivalents)",
    ]
    return "Not enough information to solve the system. Provide one of:\n  " *
           join(routes, "\n  ")
end

"Concentration unit names to the factor that converts them to mol/kg."
const UNIT_SCALES = Dict(
    "mol" => 1.0,
    "mmol" => 1.0e3,
    "umol" => 1.0e6,
    "nmol" => 1.0e9,
    "pmol" => 1.0e12,
    "fmol" => 1.0e15,
)

"""
Reduce a unit to the name `UNIT_SCALES` is keyed on, or throw listing the names it accepts.

`:umol`, `"UMOL"` and `"µmol"` all reduce to `"umol"`: a symbol or a string, any case, either
micro sign — U+00B5 from a keyboard or U+03BC from a Greek layout, neither of them the ASCII
`u` the table uses.
"""
function _normalise_unit(unit)
    name = replace(lowercase(string(unit)), 'µ' => 'u', 'μ' => 'u')
    haskey(UNIT_SCALES, name) || throw(ArgumentError("unknown unit $(repr(unit)); " *
        "expected one of " * join(sort(collect(keys(UNIT_SCALES))), ", ")))
    return name
end

# Every path normalises first, so a `KeyError` here is a bug rather than a bad argument.
_unit_multiplier(unit) = UNIT_SCALES[unit]

"A negative concentration is not a small one; it is a broken input, and NaN says so."
_clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

"Convert a concentration from the caller's unit to mol/kg."
_to_mol(v, m) = isnothing(v) ? nothing : v / m

"Convert a gas from µatm to atm."
_to_atm(v) = isnothing(v) ? nothing : v / 1e6

"""
    _to_total_scale(pHtot, pHsws, pHfree, pHNBS, env)

Reduce whichever pH the caller gave to a single total-scale value.
"""
function _to_total_scale(pHtot, pHsws, pHfree, pHNBS, env)
    if isnothing(pHtot)
        if !isnothing(pHsws)
            pHtot = pHsws - log10(env.ts_fac)
        elseif !isnothing(pHfree)
            pHtot = pHfree - log10(env.tf_fac)
        elseif !isnothing(pHNBS)
            pHtot = pHNBS + log10(env.fH) - log10(env.ts_fac)
        end
    end

    return pHtot
end

"The other three pH scales, from the solved total-scale value."
function _from_total_scale(pHtot, env)
    pHsws = pHtot + log10(env.ts_fac)
    return (pHfree = pHtot + log10(env.tf_fac),
            pHsws = pHsws,
            pHNBS = pHsws - log10(env.fH))
end

"""
Resolve CO₂, fCO₂ and pCO₂ from whichever of the three was supplied.

`fCO₂` is the one the equilibrium constant `K0` relates to CO₂; `pCO₂` reaches it through
the virial correction in `pCO₂_to_fCO₂`.
"""
function _resolve_gases(ps)
    CO₂, fCO₂, pCO₂ = ps.CO₂, ps.fCO₂, ps.pCO₂

    if isnothing(CO₂)
        if !isnothing(fCO₂)
            CO₂ = fCO₂ * ps.Ks.K0
            pCO₂ = isnothing(pCO₂) ? fCO₂_to_pCO₂(fCO₂, ps.temp_c) : pCO₂
        elseif !isnothing(pCO₂)
            fCO₂ = pCO₂_to_fCO₂(pCO₂, ps.temp_c)
            CO₂ = fCO₂ * ps.Ks.K0
        end
    end

    return (CO₂ = CO₂, fCO₂ = fCO₂, pCO₂ = pCO₂)
end

"""
Convert a supplied saturation state into the [CO₃²⁻] that produces it.

Ω is a statement about carbonate ion concentration, so this lets ΩA or ΩC stand in as one of
the two carbonate parameters — `whole_system(TA=2300.0, ΩC=5.0)` and the `carbon_system`
equivalent both work. Runs before the solvers, and for every scope that includes `:carbon`.
"""
function _omega_to_CO₃(ps, ΩA, ΩC)
    calcium = ps.Ca * ps.sal / 35.0
    if !isnothing(ΩA)
        return (CO₃ = ΩA * ps.Ks.KspA / calcium,)
    elseif !isnothing(ΩC)
        return (CO₃ = ΩC * ps.Ks.KspC / calcium,)
    end
    return NamedTuple()
end

"Saturation states of aragonite and calcite for the solved [CO₃²⁻]."
function _saturation_states(ps)
    calcium = ps.Ca * ps.sal / 35.0
    return (ΩA = ps.CO₃ * calcium / ps.Ks.KspA,
            ΩC = ps.CO₃ * calcium / ps.Ks.KspC)
end

"""
    _rescale_to_unit(ps, m, keys)

Convert every concentration in `keys` from the internal mol/kg back to the caller's unit.
"""
function _rescale_to_unit(ps, m, keys::NTuple{N,Symbol}) where {N}
    m == 1 && return ps
    scaled = map(k -> (v = getfield(ps, k); isnothing(v) ? nothing : v * m), keys)
    return merge(ps, NamedTuple{keys}(scaled))
end

"Gases are held internally in atm and reported in µatm."
_rescale_gases(ps) = (pCO₂ = ps.pCO₂ * 1e6, fCO₂ = ps.fCO₂ * 1e6)

"""
Lift the two non-numeric entries out of the computed state.

`val` on a result holds numbers only, so `Ks` and `unit` are returned separately for the
wrapper to store as their own fields. Everything else in `ps` is a concentration, a
condition or a derived quantity.

Returns `(values, Ks, unit)`.
"""
_split_metadata(ps) = (Base.structdiff(ps, NamedTuple{(:Ks, :unit)}), ps.Ks, ps.unit)



# --- One core, driven by scope -----------------------------------------------------------
#
# One procedure serves every scope. Which fields the state carries and which calculators run
# both follow from the scope, so scope is a *type* parameter: the branches on it fold away
# during specialisation, giving each scope its own compiled path from a single source.

"The state a calculation starts from: inputs converted to internal units, plus the water."
function _initial_state(scope, inputs, env, m, pHtot)
    shared = (;
        BT = env.BT, ST = env.ST, FT = env.FT, Mg = env.Mg, Ca = env.Ca, fH = env.fH,
        temp_c = inputs.temp_c, pres_bar = inputs.pres_bar, sal = inputs.sal,
        pHtot = pHtot,
        # Cleared so the inner solvers cannot read a stale scale; refilled on the way out.
        pHfree = nothing, pHsws = nothing, pHNBS = nothing,
        unit = inputs.unit, Ks = env.Ks,
    )

    carbon = :carbon in scope ? (;
        DIC = _clean(_to_mol(inputs.DIC, m)),
        TA = _to_mol(inputs.TA, m),
        CO₂ = _clean(_to_mol(inputs.CO₂, m)),
        HCO₃ = _clean(_to_mol(inputs.HCO₃, m)),
        CO₃ = _clean(_to_mol(inputs.CO₃, m)),
        PT = _clean(_to_mol(inputs.PT, m)),
        SiT = _clean(_to_mol(inputs.SiT, m)),
        H2ST = _clean(_to_mol(inputs.H2ST, m)),
        NH4T = _clean(_to_mol(inputs.NH4T, m)),
        pCO₂ = _clean(_to_atm(inputs.pCO₂)),
        fCO₂ = _clean(_to_atm(inputs.fCO₂)),
    ) : NamedTuple()

    boron = :boron in scope ? (;
        BOH₃ = _clean(_to_mol(inputs.BOH₃, m)),
        BOH₄ = _clean(_to_mol(inputs.BOH₄, m)),
    ) : NamedTuple()

    isotopes = :isotopes in scope ? (;
        δBT = inputs.δBT, δBOH₃ = inputs.δBOH₃, δBOH₄ = inputs.δBOH₄,
        ABT = inputs.ABT, ABOH₃ = inputs.ABOH₃, ABOH₄ = inputs.ABOH₄,
        alphaB = inputs.alphaB,
    ) : NamedTuple()

    return merge(shared, carbon, boron, isotopes)
end

"""
Fill in the isotope quantities the solvers work in.

δ is what people measure and A is what the equations use, so anything given as δ is converted
here, and `αB` and `δBT` fall back to modern seawater when not supplied.

`alphaB` is returned, not just used: it has to reach the solvers and the result, or the
fractionation actually applied is not the one reported.
"""
function _isotope_inputs(ps)
    δBT = (isnothing(ps.δBT) && isnothing(ps.ABT)) ? Isotopes.get_δBT() : ps.δBT
    return (;
        alphaB = something(ps.alphaB, Isotopes.get_alphaB()),
        δBT = δBT,
        ABT = !isnothing(δBT) ? Isotopes.δ11_to_A11(δBT) : ps.ABT,
        ABOH₃ = !isnothing(ps.δBOH₃) ? Isotopes.δ11_to_A11(ps.δBOH₃) : ps.ABOH₃,
        ABOH₄ = !isnothing(ps.δBOH₄) ? Isotopes.δ11_to_A11(ps.δBOH₄) : ps.ABOH₄,
    )
end

"Report the isotopes in δ notation as well as fractional abundance."
_delta_values(ps, inputs) = (;
    δBT = !isnothing(inputs.δBT) ? inputs.δBT : A11_to_δ11(ps.ABT),
    δBOH₃ = !isnothing(inputs.δBOH₃) ? inputs.δBOH₃ : A11_to_δ11(ps.ABOH₃),
    δBOH₄ = !isnothing(inputs.δBOH₄) ? inputs.δBOH₄ : A11_to_δ11(ps.ABOH₄),
)

"""
Run the speciation solvers this scope covers, in an order that works.

All of them are always run for a given scope; only the order varies, because pH is the
unknown they share — whichever subsystem the inputs determine goes first and yields pH, and
the rest follow from it. See `_first_solvable`.
"""
function _run_solvers(scope, ps)
    # Decide which solvers are needed based on the scope. This creates
    # temporary functions that have all the solver machinery *if* the
    # scope requires it, otherwise does nothing and returns the input unchanged.
    carbon(p) = :carbon in scope ? merge(p, C_calculator(; p...)) : p
    boron(p) = :boron in scope ? merge(p, B_calculator(; p...)) : p
    isotopes(p) = :isotopes in scope ? merge(p, calc_B_isotopes(; p...)) : p

    # Run the solvers in an order that works for the inputs. The first solver to run is the one
    # whose inputs determine pH; the others follow from it. Because of the block above, any solver that
    # is out of scope just returns the input unchanged, allowing the calls to be made in a chain.
    solve_first = _first_solvable(ps)
    solve_first === :boron && return isotopes(carbon(boron(ps)))
    solve_first === :isotopes && return carbon(boron(isotopes(ps)))
    solve_first === :carbon && return isotopes(boron(carbon(ps)))
    throw(ArgumentError(_underdetermined_message()))
end

"Concentrations to convert back to the caller's unit, for this scope."
_concentrations(scope) = (SHARED_CONCENTRATIONS...,
                          (:carbon in scope ? CARBON_CONCENTRATIONS : ())...,
                          (:boron in scope ? BORON_CONCENTRATIONS : ())...)

"""
    _solve_core(scope, inputs)

Solve the subsystems named by `scope`.

`inputs` carries every parameter that scope accepts, already defaulted, so this is the whole
calculation from sanitised input to finished state.
"""
function _solve_core(scope, inputs)
    m = _unit_multiplier(inputs.unit)

    # BT is converted with the other concentrations rather than inside the environment: it is
    # input sanitisation, and it matches calculate_constants's own convention, whose
    # BT_method fallbacks all return mol/kg.
    BT = _clean(_to_mol(inputs.BT, m))

    # calculate_constants returns the complete environment — constants, totals and the
    # pH-scale factors together — so a supplied `Ks` is used as-is and nothing is recombined.
    env = isnothing(inputs.Ks) ?
        calculate_constants(; inputs.temp_c, inputs.sal, inputs.pres_bar, inputs.ST,
                              inputs.FT, BT, inputs.Ca, inputs.Mg, inputs.K_method,
                              inputs.KSO4_method, inputs.BT_method, inputs.KF_method,
                              inputs.KNH3_method, inputs.Ca_method, inputs.MyAMI_mode) :
        inputs.Ks

    # Convert whichever pH the caller gave to a single total-scale value, which is what the
    # solvers work in.
    pHtot = _to_total_scale(inputs.pHtot, inputs.pHsws, inputs.pHfree, inputs.pHNBS, env)

    # Build the initial state, which is the inputs converted to internal units plus the water.
    ps = _initial_state(scope, inputs, env, m, pHtot)
    :isotopes in scope && (ps = merge(ps, _isotope_inputs(ps)))

    # If the caller gave ΩA or ΩC instead, convert it to CO₃²⁻; if the caller gave CO₂, fCO₂ or pCO₂, resolve the other two.
    if :carbon in scope
        ps = merge(ps, _resolve_gases(ps))
        ps = merge(ps, _omega_to_CO₃(ps, inputs.ΩA, inputs.ΩC))
    end

    # run the appropriate solvers in an order that works depending on the inputs and scope.
    ps = _run_solvers(scope, ps)

    # Convert the solved total-scale pH back to the other three scales.
    ps = merge(ps, _from_total_scale(ps.pHtot, env))
    
    # if isotopes are in scope, calculate the δ values from the solved fractional abundances.
    :isotopes in scope && (ps = merge(ps, _delta_values(ps, inputs)))

    # Calculate the saturation states from the solved [CO₃²⁻] and the water's Ca²⁺.
    if :carbon in scope
        ps = merge(ps, _saturation_states(ps))
    end

    # Every scope converts its concentrations back, or a boron-only result would report BT in
    # mol/kg having been given it in µmol/kg.
    ps = _rescale_to_unit(ps, m, _concentrations(scope))
    :carbon in scope && (ps = merge(ps, _rescale_gases(ps)))

    return ps
end
