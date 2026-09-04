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
    T_out = "at_collection_conditions(result; temp_c=...)",
    S_out = "nothing - a sample's salinity does not change between collection and measurement",
    P_out = "at_collection_conditions(result; pres_bar=...)",
    # Julia splatting covers this natively and composes better.
    pdict = "splatting: carbon_system(; TA=2300.0, DIC=2000.0, your_parameters...)",
    # Selected one K_method per array versus one per sample. The public entry points are
    # scalar-only, so many samples are handled by broadcasting a solver instead.
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
        if haskey(PARAMETER_DEFAULTS, name)
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

These are the union of the methods' stated domains, so being outside them means no method
covers the input — not merely that the selected one does not.
"""
const CONDITION_RANGES = (temp_c = (-5.0, 50.0), sal = (0.0, 50.0), pres_bar = (0.0, 1200.0))

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


"Warn if TA and CO₃²⁻ are supplied together as the solving pair."
function _check_ta_co3_pair(params)
    has_ta = if params isa Tuple
        :TA in params
    else
        haskey(params, :TA) && !isnothing(params.TA)
    end

    has_co3 = if params isa Tuple
        :CO3 in params || :CO₃ in params
    else
        (haskey(params, :CO3) && !isnothing(params.CO3)) ||
        (haskey(params, :CO₃) && !isnothing(params.CO₃))
    end

    if has_ta && has_co3
        @warn """
        The (TA, CO₃²⁻) parameter pair can be mathematically non-unique in certain oceanographic domains.
        Depending on your values, the numerical root-finder may fail to converge or yield non-physical solutions.
        Consult the TA vs. CO₃²⁻ stability heatmap in `examples.jl` to check if your inputs fall within the stable convergence domain.
        Still computing.
        """ maxlog=1
    end
end