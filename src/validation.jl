# Input validation for the public entry points.
#
# Julia routes any keyword a function does not declare into its `kwargs...` sink, silently.
# Every entry point in this package has such a sink, so without a guard
# `carbon_system(TA=2300, DIC=2000, tempc=2.0)` computes at the default 25 °C and returns a
# plausible wrong number. These checks reject unrecognised names, contradictory ones, and
# input that does not describe exactly one system.

"""
Collapse a keyword name to the form that survives the mistakes people actually make —
dropped underscores and inconsistent case. `tempc`, `TempC` and `temp_c` all reduce to
`tempc`, which is enough to catch the common typo without pulling in an edit-distance
implementation for the sake of it.
"""
_normalise_keyword(name::Symbol) = lowercase(replace(String(name), "_" => ""))

"Parameter names whose normalised form matches `name`, for a 'did you mean' hint."
function _similar_keywords(name::Symbol)
    target = _normalise_keyword(name)
    return [k for k in keys(PARAMETER_DEFAULTS) if _normalise_keyword(k) == target]
end

"""
Arguments this package does not accept, mapped to what to use instead.

Read only to build an error message, so a name here costs nothing until someone uses it. Each
entry is phrased as the replacement, because that is what the caller needs to type next.
"""
const RETIRED_ARGUMENTS = (
    T_in  = "temp_c",
    S_in  = "sal",
    P_in  = "pres_bar",
    T_out = "adjust(result; temp_c=...)",
    S_out = "nothing - adjust re-solves the same water, and salinity does not change with " *
            "temperature or pressure",
    P_out = "adjust(result; pres_bar=...)",
    # Julia splatting covers this natively and composes better.
    pdict = "splatting: carbon_system(; TA=2300.0, DIC=2000.0, your_parameters...)",
    # Selected one parameterisation per array versus one per sample. The public entry points
    # are scalar-only, so many samples are handled by broadcasting a solver instead.
    K_mode = "CarbonateSystem(:carbon; varying=(:TA, :DIC)), which builds a scalar " *
             "solver you can broadcast",
    # A generic pH plus `scale` says what the scale-specific arguments already say.
    pH = "pHtot, pHsws, pHfree or pHNBS - name the scale directly",
    scale = "nothing - name the scale in the argument, e.g. pHsws=8.0",
)

"""
Reject arguments this package no longer accepts.

Kept separate from [`_reject_unknown_parameters`](@ref) so a retired name gets a message naming
its replacement, rather than the 'did you mean' guess an unrecognised name gets.
"""
function _reject_retired_arguments(kwargs)
    found = [k for k in keys(kwargs) if haskey(RETIRED_ARGUMENTS, k)]
    isempty(found) && return nothing

    lines = ["$k is no longer accepted; use $(RETIRED_ARGUMENTS[k])" for k in found]

    # The condition rename needs a word of explanation that the others do not.
    conditions = [k for k in found if k in (:T_in, :S_in, :P_in, :T_out, :S_out, :P_out)]
    note = isempty(conditions) ? "" :
        "\nConditions are now named temp_c, sal and pres_bar, and a single call describes " *
        "one set of conditions."

    throw(ArgumentError("retired argument(s):\n  " * join(lines, "\n  ") * note))
end

"""
Reject parameter names this scope does not accept, naming why.

Uncertainties go through the same check, because `varying_errors` names parameters: a σ for
something that is not a parameter is a typo, and should fail when the solver is built rather
than on the first row.

`errors` gets no exemption. Exempting it makes `CarbonateSystem(:carbon; errors = (TA = 2.0,))`
legal, stored as a *parameter*, and propagating nothing — uncertainties asked for and silently
not delivered. Use `varying_errors` instead.
"""
function _reject_unknown_parameters(names, accepted, scope::Tuple{Vararg{Symbol}})
    # Allocation-free on the success path; message detail is built only when it is needed.
    all(n -> n in accepted, names) && return nothing

    lines = map(collect(n for n in names if !(n in accepted))) do name
        if name in LEGACY_METHOD_KEYWORDS
            # Reached only when the old keyword was named in `varying` or `varying_errors`;
            # given as a value, it is translated before this check.
            "$name is now a field of K_settings, which cannot vary between samples in a " *
            "solver; build one solver per setting"
        elseif haskey(PARAMETER_DEFAULTS, name)
            "$name needs the :$(getproperty(PARAMETER_SCOPE, name)) system, " *
            "but this solver's scope is $scope"
        else
            near = _similar_keywords(name)
            isempty(near) ? "$name is not a parameter" :
                "$name (did you mean $(join(near, " or "))?)"
        end
    end
    throw(ArgumentError("unrecognised parameter(s):\n  " * join(lines, "\n  ")))
end


# --- Determinacy -----------------------------------------------------------------------

# `PARAMETER_GROUPS` lives in core.jl, beside the solve-order dispatch that reads the same
# table. Counting parameters rather than groups would treat `pHtot=8.1, pHsws=8.0` as a
# solvable system, and it is not — it is a contradiction.
#
# The boron entries matter because pH is not only a carbon observable: `BT` with one of
# `BOH₃`/`BOH₄` fixes it, as does `δBT` with one of `δBOH₃`/`δBOH₄`. Only `whole_system`
# accepts them; `carbon_system` takes `BT` alone, so for that function the carbon groups are
# the whole story and any boron speciation name is rejected as unrecognised first.

"How many of a group's member names `inputs` supplies."
function _supplied_in_group(inputs::NamedTuple, members)
    n = 0
    for name in members
        haskey(inputs, name) && !isnothing(inputs[name]) && (n += 1)
    end
    return n
end

"""
The member names of `members` that `inputs` supplies.

Only called to build an error message, never on the path of a successful calculation, so it
is free to allocate.
"""
_supplied_names(inputs::NamedTuple, members) =
    [name for name in members if haskey(inputs, name) && !isnothing(inputs[name])]

"The one name supplied for each group that has exactly one, for error messages."
function _supplied_names(inputs::NamedTuple)
    names = Symbol[]
    for (_, members) in pairs(PARAMETER_GROUPS)
        append!(names, _supplied_names(inputs, members))
    end
    return names
end

# Scope stands in for the entry point in validation messages.
_entry_name(scope::Tuple{Vararg{Symbol}}) =
    scope === (:carbon,) ? :carbon_system :
    scope === (:carbon, :boron, :isotopes) ? :whole_system : Symbol("CarbonateSystem", scope)

"""
    _check_determinacy(inputs, scope; require_two)

Reject input that does not describe exactly one system.

Two independent constraints determine the carbonate system; fewer leaves it unsolvable, more
over-determines it, and both have to be errors. Over-determined input is the more dangerous
of the two: `carbon_system(TA=2300, DIC=2000, pHtot=7.0)` would otherwise discard one of the
three and recompute it — returning `DIC = 2445.76` with nothing in the result to show that a
supplied measurement was dropped.

`require_two` is false for `whole_system`, whose scope can be determined by boron or isotope
parameters instead. Its core reports the under-determined case, listing every route to a
solution; duplicating that here would mean two places to keep in step.
"""
function _check_determinacy(inputs::NamedTuple, scope::Tuple{Vararg{Symbol}};
                            require_two::Bool)
    # Counted without allocating, because this runs on every call and answers a yes/no
    # question. Building the collections needed to *describe* the problem allocates, so that
    # work is pushed onto the error paths, where it happens once.
    groups = 0
    for (group, members) in pairs(PARAMETER_GROUPS)
        supplied = _supplied_in_group(inputs, members)
        supplied == 0 && continue

        if supplied > 1
            present = _supplied_names(inputs, members)
            throw(ArgumentError(
                "$(_entry_name(scope)) was given $(join(present, ", ")), which all " *
                "describe $group. Supply one of them.\nThese are the same quantity " *
                "expressed differently rather than independent measurements, so all but " *
                "one would be silently ignored."
            ))
        end

        groups += 1
    end

    if groups > 2
        throw(ArgumentError(
            "$(_entry_name(scope)) was given $groups parameters " *
            "($(join(_supplied_names(inputs), ", "))); two determine the system.\n" *
            "Drop one: the extra value would be discarded and recomputed, with nothing in " *
            "the result to show that the supplied measurement disagreed."
        ))
    end

    if require_two && groups < 2
        names = _supplied_names(inputs)
        throw(ArgumentError(
            "$(_entry_name(scope)) needs two parameters to solve the system, but was " *
            "given $(isempty(names) ? "none" : "only $(only(names))").\nSupply two of: " *
            "TA, DIC, pH (or pHtot/pHsws/pHfree/pHNBS), CO₂ (or pCO₂/fCO₂), HCO₃, " *
            "CO₃ (or ΩA/ΩC)."
        ))
    end

    return nothing
end


# --- Conditions ------------------------------------------------------------------------

"""
Ranges outside which every available parameterisation is an extrapolation rather than a fit.

Temperature and salinity are the union of the K1/K2 parameterisations' fitted ranges, built
from the registry ([`K1_K2_RANGES`](@ref)), so being outside them means no method covers the
input — not merely that the selected one does not; `_check_K_settings` covers that.
"""
const CONDITION_RANGES = (
    temp_c = (minimum(r.temp_c[1] for r in values(K1_K2_RANGES)),
              maximum(r.temp_c[2] for r in values(K1_K2_RANGES))),
    sal = (minimum(r.sal[1] for r in values(K1_K2_RANGES)),
           maximum(r.sal[2] for r in values(K1_K2_RANGES))),
    pres_bar = (0.0, 1200.0),
)

"""
    _check_conditions(temp_c, sal, pres_bar)

Reject impossible conditions, and warn about implausible ones.

Impossible means the chemistry is undefined: a negative salinity makes the ionic-strength
terms in the constants meaningless, and would otherwise surface as a `DomainError` from inside
`sqrt` with nothing naming the input responsible.

Out-of-range but physical values only warn. The package should still compute for a brine or a
hydrothermal vent, where refusing would be worse than extrapolating. `maxlog=1` keeps a large
batch readable, since this runs on every call.
"""
function _check_conditions(temp_c, sal, pres_bar)
    # Arrays and ForwardDiff.Dual values reach here on some paths; only plain numbers can be
    # meaningfully compared, and the public API is scalar-only regardless.
    all(x -> x isa Real, (temp_c, sal, pres_bar)) || return nothing

    for (name, value) in ((:temp_c, temp_c), (:sal, sal), (:pres_bar, pres_bar))
        isfinite(value) || throw(ArgumentError("$name must be finite, got $value"))
    end

    sal < 0 && throw(ArgumentError(
        "sal must not be negative, got $sal.\nSalinity below zero makes the ionic-strength " *
        "terms in the equilibrium constants undefined, which surfaces as a DomainError " *
        "from inside sqrt rather than as anything naming the input."))
    pres_bar < 0 && throw(ArgumentError("pres_bar must not be negative, got $pres_bar"))

    for (name, value) in ((:temp_c, temp_c), (:sal, sal), (:pres_bar, pres_bar))
        low, high = getproperty(CONDITION_RANGES, name)
        if value < low || value > high
            @warn "$name = $value is outside $low to $high, where every available " *
                  "parameterisation is an extrapolation. Still computing." maxlog = 1
        end
    end

    return nothing
end


# --- Suitability of the chosen constants --------------------------------------------------

"Units for the conditions in the suitability warnings."
_condition_unit(name::Symbol) = name === :temp_c ? " °C" : ""

"""
    _check_K_settings(settings, temp_c, sal, Ca, Mg)

Warn when the chosen K1/K2 parameterisation does not suit this water: a temperature or
salinity outside the range it was fitted over, or a non-modern seawater composition that only
KGen corrects for.

Warns rather than refusing, like `_check_conditions`: an extrapolation is often the best
available, and the caller is the one who can judge that. The choice itself is never changed —
picking a parameterisation from the conditions is what `K_method = "automatic"` did, and it
made the constants step discontinuously between samples.

`maxlog=1` per method, condition and side, so a broadcast over a cruise warns once for each
way the data falls outside the range rather than once per sample. A broadcast has no batch to
count over, which is why the warning cannot say how many samples were affected.
"""
function _check_K_settings(settings::KSettings, temp_c, sal, Ca, Mg)
    # Duals reach here on the AD path, and only plain numbers compare meaningfully.
    all(x -> x isa Real, (temp_c, sal)) || return nothing
    method = settings.K1_K2_method
    fit = K1_K2_RANGES[method]

    for (name, value) in ((:temp_c, temp_c), (:sal, sal))
        low, high = getproperty(fit, name)
        low <= value <= high && continue
        side = value < low ? :below : :above
        unit = _condition_unit(name)
        fitted = low == high ? "$name = $low$unit only" : "$low to $high$unit"
        @warn "You currently have selected \"$method\" for K1_K2_method, however this " *
              "method is not recommended for $name = $value$unit: it was fitted over " *
              "$fitted ($(fit.source)). Still computing. Shown once for each method and " *
              "side of its range; the README tabulates every K1_K2_method's range." maxlog=1 _id=(
                  :K1_K2_range, method, name, side)
    end

    # Only Kgen corrects the constants for Ca and Mg; every other parameterisation assumes
    # modern seawater whatever composition the saturation states are computed with.
    method == "KGen" && return nothing
    non_modern(x, modern) = x isa Real && x != modern
    if non_modern(Ca, MODERN_CALCIUM) || non_modern(Mg, MODERN_MAGNESIUM)
        @warn "You currently have selected \"$method\" for K1_K2_method, however this " *
              "method is not recommended for a non-modern seawater composition " *
              "(Ca = $(something(Ca, MODERN_CALCIUM)), Mg = $(something(Mg, MODERN_MAGNESIUM)) " *
              "mol/kg): it assumes modern seawater, and only \"KGen\" corrects the " *
              "constants for Ca and Mg. Still computing." maxlog=1 _id=(:K1_K2_composition, method)
    end
    return nothing
end
