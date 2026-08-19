# The solver: one type, callable two ways.

"""
    CarbonateSystem{scope, varying, varying_errors}

A carbonate system solver.

`scope` names the subsystems it computes, `varying` the parameters supplied positionally,
and `varying_errors` the measurement uncertainties supplied positionally after them (both
empty for the keyword form). All three are type parameters, so the branches on them resolve
during specialisation.

Two ways to call one object:

```julia
sys = CarbonateSystem(:carbon)
sys(TA = 2300.0, DIC = 2000.0)                      # keyword

fast = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c))
fast.([2300.0, 2310.0], [2000.0, 1990.0], [25.0, 2.0])   # positional, so it broadcasts
```

[`carbon_system`](@ref), [`whole_system`](@ref), [`boron_system`](@ref) and
[`boron_isotopes`](@ref) are presets over the keyword form.
"""
struct CarbonateSystem{scope, varying, varying_errors, D<:NamedTuple, S<:NamedTuple}
    defaults::D
    settings::S
end

"""
    CarbonateSystem(scope::Symbol...; varying = (), varying_errors = (), settings...)

Build a solver.

`scope` is some combination of `:carbon`, `:boron` and `:isotopes`; `varying` names the
parameters that will be given positionally, in order; everything else is a fixed `setting`,
stated once and unable to drift between calls.

`varying_errors` names measurement uncertainties supplied the same way, positionally after
the varying parameters. Every uncertainty is varying because a scalar broadcasts against a
vector, so a σ that is the same for every sample needs no separate mechanism — just pass the
number:

```julia
solve = CarbonateSystem(:carbon; varying = (:TA, :DIC), varying_errors = (:TA, :DIC))

TA, DIC = [2300.0, 2310.0], [2000.0, 1990.0]
solve.(TA, DIC, [2.0, 3.0], [1.5, 2.0])   # a σ column per sample, in the order specified in the construcor
solve.(TA, DIC, 2.0, 2.0)                 # or one σ for all of them
```

Everything knowable from names alone is checked here rather than per call — an unrecognised
parameter, a scope that cannot use it, a set that does not determine the system, an
uncertainty named for a value the solver does not have. A broadcast over a million rows
therefore fails at construction, not on row 700,000.
"""
function CarbonateSystem(scope::Symbol...; varying = (), varying_errors = (), settings...)
    known = (:carbon, :boron, :isotopes)
    for s in scope
        s in known || throw(ArgumentError("unknown scope :$s; expected some of $known"))
    end
    isempty(scope) && throw(ArgumentError("name at least one scope, e.g. :carbon"))

    varying = Tuple(varying)
    varying_errors = Tuple(varying_errors)
    fixed = NamedTuple(settings)

    for (names, what) in ((varying, "varying"), (varying_errors, "varying_errors"))
        repeated = unique(n for n in names if count(==(n), names) > 1)
        isempty(repeated) || throw(ArgumentError(
            "$(join(repeated, ", ")) named more than once in `$what`"))
    end

    overlap = [n for n in varying if haskey(fixed, n)]
    isempty(overlap) || throw(ArgumentError(
        "$(join(overlap, ", ")) given both as a varying parameter and as a fixed setting; " *
        "it has to be one or the other"))

    # An uncertainty propagates from a measurement, so it needs a value to attach to.
    for name in varying_errors
        name in varying || haskey(fixed, name) || throw(ArgumentError(
            "$name is in `varying_errors`, but this solver has no value for it. Put $name " *
            "in `varying`, or give it as a fixed setting."))
    end

    _reject_unknown_parameters(
        (varying..., varying_errors..., keys(fixed)...),
        _accepted_parameters(scope),
        scope
        )

    # Resolved once here and stored on the solver. Building a NamedTuple from a
    # runtime-computed key set is the expensive pattern in this package — doing it per call
    # rather than per solver costs about 50% of a whole calculation.
    defaults = NamedTuple{_accepted_parameters(scope)}(PARAMETER_DEFAULTS)

    # Values are stand-ins: the determinacy check cares which parameters are present, not
    # what they hold.
    if !isempty(varying) || !isempty(fixed)
        prototype = merge(fixed, NamedTuple{varying}(ntuple(_ -> 1.0, length(varying))))
        isempty(varying) ||
            _check_determinacy(prototype, scope; require_two = scope === (:carbon,))
    end

    return CarbonateSystem{scope, varying, varying_errors, typeof(defaults),
                           typeof(fixed)}(defaults, fixed)
end

scope_of(::CarbonateSystem{scope}) where {scope} = scope
varying_of(::CarbonateSystem{scope, varying}) where {scope, varying} = varying

"""
Solve, given parameters by keyword.

There is no `errors` keyword here. Uncertainty propagation is reached by building a solver
with `varying_errors` and calling it positionally, which keeps the everyday call — the one
`carbon_system` and friends wrap — to just the chemistry.
"""
(sys::CarbonateSystem)(; kwargs...) = _run(sys, NamedTuple(kwargs), nothing)

"""
Solve, given the `varying` parameters positionally, then the `varying_errors`.

This is the form that broadcasts, because Julia will not broadcast over keyword arguments,
and the only form that propagates uncertainties.
"""
function (sys::CarbonateSystem{scope, varying, varying_errors})(
        values::Vararg{Any, N}) where {scope, varying, varying_errors, N}
    # `varying_errors` alone is legitimate: fixed water, varying uncertainty.
    isempty(varying) && isempty(varying_errors) && throw(ArgumentError(
        "this solver takes keyword parameters; build it with `varying = (:TA, :DIC)` or " *
        "`varying_errors = (:DIC,)` to call it positionally"))

    expected = length(varying) + length(varying_errors)
    N == expected || throw(ArgumentError(
        "this solver takes $expected values, in the order " *
        join((string.(varying)..., ("σ($e)" for e in varying_errors)...), ", ") *
        "; got $N"))

    parameters = NamedTuple{varying}(values[1:length(varying)])
    errors = isempty(varying_errors) ? nothing :
             NamedTuple{varying_errors}(values[length(varying) + 1:end])

    return _run(sys, parameters, errors)
end

"Validate, fill in defaults, compute and package — the one path both call forms take."
function _run(sys::CarbonateSystem{scope}, supplied_parameters, errors) where {scope}
    supplied = merge(sys.settings, supplied_parameters)
    _reject_retired_arguments(supplied)
    _reject_unknown_parameters(keys(supplied), keys(sys.defaults), scope)

    inputs = merge(sys.defaults, supplied)

    _check_determinacy(inputs, scope; require_two = scope === (:carbon,))
    _check_conditions(inputs.temp_c, inputs.sal, inputs.pres_bar)

    return _solve(scope, inputs, errors)
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

# Built once at load, not inside the presets above. Constructing a solver resolves its
# accepted parameter set and defaults; doing that per call costs about 50% of a calculation.
const CARBON_SYSTEM = CarbonateSystem(:carbon)
const WHOLE_SYSTEM = CarbonateSystem(:carbon, :boron, :isotopes)
const BORON_SYSTEM = CarbonateSystem(:boron, :isotopes)
const BORON_ISOTOPES = CarbonateSystem(:isotopes)

# There is deliberately no `carbon_solver` / `whole_solver` pair: a solver is just
# `CarbonateSystem(scope; varying = …)`, which reads the same for every scope. Naming two of
# the eight scope/varying combinations would imply "solver" is a separate concept.

function Base.show(io::IO, sys::CarbonateSystem{scope, varying}) where {scope, varying}
    print(io, "CarbonateSystem(", join((":$s" for s in scope), ", "), ")")
    isempty(varying) || print(io, " varying ", join((":$v" for v in varying), ", "))
    isempty(sys.settings) ||
        print(io, " with ", join(("$k=$(repr(v))" for (k, v) in pairs(sys.settings)), ", "))
end
