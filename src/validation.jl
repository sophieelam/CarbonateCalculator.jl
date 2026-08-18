# Input validation for the public entry points.
#
# Julia routes any keyword a function does not declare into its `kwargs...` sink, silently.
# Every entry point in this package has such a sink, so `carbon_system(TA=2300, DIC=2000,
# tempc=2.0)` computes at the default 25 °C and returns a plausible wrong number. That is
# the same trap `_reject_retired_arguments` was added for during the conditions refactor;
# this file generalises the guard to every unrecognised name.

"""
    _declared_keywords(entry_point) -> Vector{Symbol}

The keyword names `entry_point` actually declares, read off its method signature.

Used only to suggest a correction once a call has already failed, so it is never on the hot
path, and reading the signature rather than maintaining a parallel list means the two cannot
drift apart. The `kwargs...` sink itself appears as a name ending in `...` and is dropped.
"""
function _declared_keywords(entry_point)
    names = Symbol[]
    for method in methods(entry_point)
        append!(names, Base.kwarg_decl(method))
    end
    filter!(name -> !endswith(String(name), "..."), names)
    return sort!(unique!(names))
end

"""
Collapse a keyword name to the form that survives the mistakes people actually make —
dropped underscores and inconsistent case. `tempc`, `TempC` and `temp_c` all reduce to
`tempc`, which is enough to catch the common typo without pulling in an edit-distance
implementation for the sake of it.
"""
_normalise_keyword(name::Symbol) = lowercase(replace(String(name), "_" => ""))

"""
Parameter names whose normalised form matches `name`, for a 'did you mean' hint.

Reads `PARAMETER_DEFAULTS` when the entry point is a solver, and the method signature when it
is a plain function, so both paths get the same suggestion.
"""
function _similar_keywords(name::Symbol, entry_point)
    target = _normalise_keyword(name)
    candidates = entry_point isa Tuple ? collect(keys(PARAMETER_DEFAULTS)) :
                                         _declared_keywords(entry_point)
    return [k for k in candidates if _normalise_keyword(k) == target]
end

"""
    _reject_unknown_arguments(kwargs, entry_point)

Reject keywords the calculation does not understand.

`kwargs` is the entry point's keyword sink. Julia routes only *unmatched* keywords there, so
anything present is by definition a name `entry_point` does not accept — the check needs no
list of valid names, and cannot fall out of step with the signature.

Retired condition names are handled first so they keep their more specific message, which
names the replacement.
"""
function _reject_unknown_arguments(kwargs, entry_point)
    _reject_retired_arguments(kwargs)
    isempty(kwargs) && return nothing

    lines = map(collect(keys(kwargs))) do name
        suggestions = _similar_keywords(name, entry_point)
        isempty(suggestions) ? "$name" : "$name (did you mean $(join(suggestions, " or "))?)"
    end

    throw(ArgumentError(
        "unrecognised argument(s) to $(_entry_name(entry_point)):\n  " * join(lines, "\n  ") *
        "\nUnrecognised keywords used to be absorbed silently, which meant a typo returned " *
        "a plausible result computed at the defaults."
    ))
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

"""
    _check_determinacy(inputs, entry_point; require_two)

Reject input that does not describe exactly one system.

Two independent constraints determine the carbonate system; fewer leaves it unsolvable, more
over-determines it. Both used to fail silently. An under-determined call surfaced as
`UndefVarError: H not defined in local scope`, naming an internal variable rather than the
mistake. An over-determined one quietly discarded a supplied value and recomputed it, so
`carbon_system(TA=2300, DIC=2000, pHtot=7.0)` returned `DIC = 2445.76` without a word — the
worse of the two, because nothing in the result records that a measurement was dropped.

`require_two` is false for `whole_system`: its core already reports the under-determined
case usefully, listing the boron and isotope routes to a solution, and duplicating that
logic here would mean two places to keep in step.
"""
function _check_determinacy(inputs::NamedTuple, entry_point; require_two::Bool)
    # Counted without allocating, because this runs on every call. The first version built a
    # Vector{Pair{Symbol,Vector{Symbol}}} to hold the same information and cost 7.6 µs - 24%
    # of a whole `carbon_system` call, for a check that answers a yes/no question. Message
    # construction is pushed onto the error paths, where allocation costs nothing.
    groups = 0
    for (group, members) in pairs(PARAMETER_GROUPS)
        supplied = _supplied_in_group(inputs, members)
        supplied == 0 && continue

        if supplied > 1
            present = _supplied_names(inputs, members)
            throw(ArgumentError(
                "$(_entry_name(entry_point)) was given $(join(present, ", ")), which all " *
                "describe $group. Supply one of them.\nThese are the same quantity " *
                "expressed differently rather than independent measurements, so all but " *
                "one would be silently ignored."
            ))
        end

        groups += 1
    end

    if groups > 2
        throw(ArgumentError(
            "$(_entry_name(entry_point)) was given $groups parameters " *
            "($(join(_supplied_names(inputs), ", "))); two determine the system.\n" *
            "Drop one: the extra value would be discarded and recomputed, with nothing in " *
            "the result to show that the supplied measurement disagreed."
        ))
    end

    if require_two && groups < 2
        names = _supplied_names(inputs)
        throw(ArgumentError(
            "$(_entry_name(entry_point)) needs two parameters to solve the system, but was " *
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

Negative salinity used to surface as a `DomainError` from a `sqrt` several frames down, with
nothing naming the input responsible. Out-of-range but physical values only warn: the
package should still compute for a brine or a hydrothermal vent, and refusing would be worse
than extrapolating. `maxlog=1` because this sits on the path of every call, and a warning
per row would make a large batch unreadable.
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
