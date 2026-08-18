# The solver: one type, callable two ways.

"""
    CarbonateSystem{scope, varying}

A carbonate system solver.

`scope` names the subsystems it computes, and `varying` the parameters supplied positionally
(empty for the keyword form). Both are type parameters, so the branches on them resolve
during specialisation.

Two ways to call one object:

```julia
sys = CarbonateSystem(:carbon)
sys(TA = 2300.0, DIC = 2000.0)                      # keyword

fast = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c))
fast.(df.talk, df.tco2, df.temperature)             # positional, so it broadcasts
```

[`carbon_system`](@ref) and [`whole_system`](@ref) are presets over the first form;
[`carbon_solver`](@ref) and [`whole_solver`](@ref) over the second.
"""
struct CarbonateSystem{scope, varying, D<:NamedTuple, S<:NamedTuple}
    defaults::D
    settings::S
end

"""
    CarbonateSystem(scope::Symbol...; varying = (), settings...)

Build a solver.

`scope` is some combination of `:carbon`, `:boron` and `:isotopes`; `varying` names the
parameters that will be given positionally, in order; everything else is a fixed `setting`,
stated once and unable to drift between calls.

Everything knowable from names alone is checked here rather than per call — an unrecognised
parameter, a scope that cannot use it, a set that does not determine the system. A broadcast
over a million rows therefore fails at construction, not on row 700,000.
"""
function CarbonateSystem(scope::Symbol...; varying = (), settings...)
    known = (:carbon, :boron, :isotopes)
    for s in scope
        s in known || throw(ArgumentError("unknown scope :$s; expected some of $known"))
    end
    isempty(scope) && throw(ArgumentError("name at least one scope, e.g. :carbon"))

    varying = Tuple(varying)
    fixed = NamedTuple(settings)

    duplicates = unique(n for n in varying if count(==(n), varying) > 1)
    isempty(duplicates) || throw(ArgumentError(
        "$(join(duplicates, ", ")) named more than once in `varying`"))

    overlap = [n for n in varying if haskey(fixed, n)]
    isempty(overlap) || throw(ArgumentError(
        "$(join(overlap, ", ")) given both as a varying parameter and as a fixed setting; " *
        "it has to be one or the other"))

    _reject_unknown_parameters(
        (varying..., keys(fixed)...), 
        _accepted_parameters(scope),
        scope
        )

    # Resolved once, here, and stored. Computing the accepted key set per call is what made
    # an earlier version of this 50% slower than the wrappers it replaced: it rebuilt a
    # NamedTuple from a runtime-computed key set on every call, the exact pattern phase 6c
    # removed from `_settings`, `_rescale_to_unit` and `_finalise`.
    defaults = NamedTuple{_accepted_parameters(scope)}(PARAMETER_DEFAULTS)

    # Values are stand-ins: the determinacy check cares which parameters are present, not
    # what they hold.
    if !isempty(varying) || !isempty(fixed)
        prototype = merge(fixed, NamedTuple{varying}(ntuple(_ -> 1.0, length(varying))))
        isempty(varying) ||
            _check_determinacy(prototype, scope; require_two = scope === (:carbon,))
    end

    return CarbonateSystem{scope, varying, typeof(defaults), typeof(fixed)}(defaults, fixed)
end

scope_of(::CarbonateSystem{scope}) where {scope} = scope
varying_of(::CarbonateSystem{scope, varying}) where {scope, varying} = varying

"""
Solve, given parameters by keyword.

The whole user-facing path in one place, for every scope: validate, fill in defaults,
compute, package.
"""
function (sys::CarbonateSystem{scope})(; errors = nothing, kwargs...) where {scope}
    supplied = merge(sys.settings, NamedTuple(kwargs))
    _reject_retired_arguments(supplied)
    _reject_unknown_parameters(keys(supplied), keys(sys.defaults), scope)

    inputs = merge(sys.defaults, supplied)

    _check_determinacy(inputs, scope; require_two = scope === (:carbon,))
    _check_conditions(inputs.temp_c, inputs.sal, inputs.pres_bar)
    _check_error_names(errors, inputs, scope)

    return _solve(scope, inputs, errors)
end

"""
Solve, given the `varying` parameters positionally.

This is the form that broadcasts, because Julia will not broadcast over keyword arguments.
"""
function (sys::CarbonateSystem{scope, varying})(values::Vararg{Any, N}) where {scope, varying, N}
    isempty(varying) && throw(ArgumentError(
        "this solver takes keyword parameters; build it with `varying = (:TA, :DIC)` to " *
        "call it positionally"))
    N == length(varying) || throw(ArgumentError(
        "this solver takes $(length(varying)) values, in the order " *
        "$(join(varying, ", ")); got $N"))

    return sys(; NamedTuple{varying}(values)...)
end

"Run the pipeline for `scope` and package the result."
function _solve(scope, inputs, errors)
    provided = [k for k in CARBONATE_PARAMETERS
                if haskey(inputs, k) && !isnothing(inputs[k])]

    # `propagate_errors` differentiates a keyword function, so the scope is closed over.
    core(; kwargs...) = _solve_core(scope, NamedTuple(kwargs))

    isnothing(errors) && return CarbonateResult(_solve_core(scope, inputs), nothing,
                                                provided, _settings(inputs), inputs, scope,
                                                nothing)

    result = propagate_errors(core; inputs = inputs, errors = errors)
    return CarbonateResult(result.val, result.err, provided, _settings(inputs), inputs,
                           scope, errors)
end

"Which carbonate parameters `show` treats as inputs rather than results."
const CARBONATE_PARAMETERS = (:TA, :DIC, :pHtot, :pCO₂, :fCO₂, :CO₃, :HCO₃)

"Reject parameter names this scope does not accept, naming why."
function _reject_unknown_parameters(names, accepted, scope::Tuple{Vararg{Symbol}})
    # Allocation-free on the success path; message detail is built only when it is needed.
    all(n -> n in accepted || n === :errors, names) && return nothing

    lines = map(collect(n for n in names if !(n in accepted) && n !== :errors)) do name
        if haskey(PARAMETER_DEFAULTS, name)
            "$name needs the :$(getproperty(PARAMETER_SCOPE, name)) system, " *
            "but this solver's scope is $scope"
        else
            near = [k for k in keys(PARAMETER_DEFAULTS)
                    if _normalise_keyword(k) == _normalise_keyword(name)]
            isempty(near) ? "$name is not a parameter" :
                "$name (did you mean $(join(near, " or "))?)"
        end
    end
    throw(ArgumentError("unrecognised parameter(s):\n  " * join(lines, "\n  ")))
end

# Scope stands in for the entry point in validation messages.
_entry_name(scope::Tuple{Vararg{Symbol}}) =
    scope === (:carbon,) ? :carbon_system :
    scope === (:carbon, :boron, :isotopes) ? :whole_system : Symbol("CarbonateSystem", scope)

# boron_system and boron_isotopes are still plain functions, and reach the same checks.
_entry_name(entry_point) = nameof(entry_point)


# --- The user-facing presets ------------------------------------------------------------

"""
    carbon_system(; TA, DIC, temp_c, sal, ...)

Solve the carbonate system. A preset over [`CarbonateSystem`](@ref) with scope `(:carbon,)`.

Every parameter it accepts, and every default, comes from `PARAMETER_DEFAULTS`.

```julia
carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0).pHtot
```
"""
carbon_system(; kwargs...) = CARBON_SYSTEM(; kwargs...)

"""
    whole_system(; TA, DIC, δBT, temp_c, sal, ...)

Solve the carbonate, boron and boron-isotope systems together. A preset over
[`CarbonateSystem`](@ref) with scope `(:carbon, :boron, :isotopes)`.

Accepts everything `carbon_system` does, plus the boron speciation and isotope parameters.
"""
whole_system(; kwargs...) = WHOLE_SYSTEM(; kwargs...)

"""
    boron_system(; pHtot, BT, BOH₃, BOH₄, δBT, ...)

Solve boron speciation and its isotopes, without the carbonate system. A preset over
[`CarbonateSystem`](@ref) with scope `(:boron, :isotopes)`.

Carbonate parameters are not accepted: with no carbon in scope, `DIC` is an unrecognised
parameter rather than a value quietly ignored.
"""
boron_system(; kwargs...) = BORON_SYSTEM(; kwargs...)

"""
    boron_isotopes(; pHtot, δBT, δBOH₃, δBOH₄, alphaB, ...)

Solve the boron isotope system alone. A preset over [`CarbonateSystem`](@ref) with scope
`(:isotopes,)`.
"""
boron_isotopes(; kwargs...) = BORON_ISOTOPES(; kwargs...)

# Built once, at load. Constructing a solver resolves its accepted parameter set and its
# defaults, which is cheap but not free — doing it inside the preset meant paying that on
# every call, and made these 50% slower than the hand-written wrappers they replaced.
const CARBON_SYSTEM = CarbonateSystem(:carbon)
const WHOLE_SYSTEM = CarbonateSystem(:carbon, :boron, :isotopes)
const BORON_SYSTEM = CarbonateSystem(:boron, :isotopes)
const BORON_ISOTOPES = CarbonateSystem(:isotopes)

# There is deliberately no `carbon_solver` / `whole_solver` pair. They were only
# `CarbonateSystem(scope; varying = …)` with the scope baked in, and they existed for two of
# the four presets — no `boron_solver`, no `isotopes_solver` — which is the tell that
# "solver" was not a separate concept, just two of eight arbitrary combinations. Building one
# is `CarbonateSystem(:carbon; varying = (:TA, :DIC))`, which reads the same for every scope.

function Base.show(io::IO, sys::CarbonateSystem{scope, varying}) where {scope, varying}
    print(io, "CarbonateSystem(", join((":$s" for s in scope), ", "), ")")
    isempty(varying) || print(io, " varying ", join((":$v" for v in varying), ", "))
    isempty(sys.settings) ||
        print(io, " with ", join(("$k=$(repr(v))" for (k, v) in pairs(sys.settings)), ", "))
end
