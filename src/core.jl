# The pipeline shared by `carbon_system_core` and `whole_system_core`.
#
# The stages are named functions rather than one flag-driven core because the cores marked
# them with comments (`# --- THE WAY IN ---`, `# Converting values back into their original
# units`), which is the signal to split.

# --- What determines the system -----------------------------------------------------
#
# One table, read by two things that must agree: the input validation in `validation.jl`,
# and the solve-order dispatch in `whole_system_core`. They used to be separate encodings —
# `PARAMETER_GROUPS` on one side, three hand-written `count(!isnothing, ...)` expressions
# and a hand-maintained error message on the other.

"""
Parameters grouped by the degree of freedom each constrains.

Two members of one group are not two measurements. pH on four scales is one measurement
expressed four ways; ΩA and ΩC are both statements about [CO₃²⁻]; pCO₂ and fCO₂ are both
statements about dissolved CO₂; δ and A are the same isotope value in two notations.

**`BT`, `δBT` and `ABT` are deliberately absent.** They are totals with defaults — `BT` from
salinity, `δBT` from modern seawater — so supplying one constrains nothing on its own. It is
the *speciated* member that carries the information.
"""
const PARAMETER_GROUPS = (
    pH    = (:pHtot, :pHsws, :pHfree, :pHNBS),
    CO₂   = (:CO₂, :pCO₂, :fCO₂),
    CO₃   = (:CO₃, :ΩA, :ΩC),
    HCO₃  = (:HCO₃,),
    DIC   = (:DIC,),
    TA    = (:TA,),
    BOH₃  = (:BOH₃,),
    BOH₄  = (:BOH₄,),
    δBOH₃ = (:δBOH₃, :ABOH₃),
    δBOH₄ = (:δBOH₄, :ABOH₄),
)

"""
Which subsystem each group constrains.

`pH` is `:shared` because it is the unknown that links the three: fix it in any one
subsystem and the other two follow. That is the whole basis of the solve order below.
"""
const GROUP_SUBSYSTEM = (pH = :shared, CO₂ = :carbon, CO₃ = :carbon, HCO₃ = :carbon,
                         DIC = :carbon, TA = :carbon, BOH₃ = :boron, BOH₄ = :boron,
                         δBOH₃ = :isotopes, δBOH₄ = :isotopes)

"How many independent constraints a subsystem needs before it can be solved for pH."
const CONSTRAINTS_NEEDED = (carbon = 2, boron = 1, isotopes = 1)

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
    # The common case lands here: two carbonate parameters and no boron speciation, so
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

_unit_multiplier(unit) = get(UNIT_SCALES, unit, 1.0)

"A negative concentration is not a small one; it is a broken input, and NaN says so."
_clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

"Convert a concentration from the caller's unit to mol/kg."
_to_mol(v, m) = isnothing(v) ? nothing : v / m

"Convert a gas from µatm to atm."
_to_atm(v) = isnothing(v) ? nothing : v / 1e6

"""
Every concentration held internally in mol/kg that must be converted back to the caller's
unit on the way out.
"""
const CARBON_CONCENTRATIONS = (:DIC, :TA, :BT, :CO₂, :HCO₃, :CO₃, :PT, :SiT,
                               :CAlk, :BAlk, :PAlk, :OH, :SiAlk, :HSO₄, :Hfree, :HF,
                               :Alk_H2S, :Alk_NH3)

const WHOLE_CONCENTRATIONS = (CARBON_CONCENTRATIONS..., :BOH₃, :BOH₄)

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

Ω is a statement about carbonate ion concentration, so this lets ΩA or ΩC stand in as one
of the two carbonate parameters. It lived only in the carbon core, which is why
`whole_system(TA=2300.0, ΩC=5.0)` used to fail with "Impossible! You haven't provided
enough information" while `carbon_system` solved the same input.
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
