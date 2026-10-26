# Re-solving a result at other conditions.
#
# A calculation (`carbon_system`, `whole_system`, or a `CarbonateSystem`) solves the system at
# the temperature and pressure it was given. `adjust` re-solves the same water at a different
# temperature and/or pressure, carrying its conservative totals across. Reporting a sample
# measured in the laboratory at the conditions it was collected at is one use; others include
# temperature-normalising a dataset, or asking how a water mass changes as it warms or sinks.
#
# Salinity is not a parameter: changing it changes which water this is, not its conditions.

"""
Carbonate parameters that describe the *solved state* rather than the water itself. They
must be cleared before re-solving, or the new calculation would be driven by values that only
hold at the original conditions.
"""
const DERIVED_PARAMETERS = (:pH, :pHtot, :pHsws, :pHfree, :pHNBS,
                            :CO₂, :HCO₃, :CO₃, :pCO₂, :fCO₂,
                            :ΩC, :ΩA, :Ks)

"""
The calculation that produced a result, as a keyword function.

A result records the *scope* it was solved with, so re-solving needs no guesswork about which
entry point produced it, and works for every scope.
"""
_core_for(scope::Tuple{Vararg{Symbol}}) = (; kwargs...) -> _solve_core(scope, NamedTuple(kwargs))

"""
Build the argument list that re-solves the water described by `inputs` at new conditions, with
TA and DIC as given.

TA and DIC are the conservative pair: they are what carries the chemistry across a change of
temperature or pressure. Everything else — the totals, the seawater composition, and every
setting — is inherited from `inputs`, which keeps the two calculations consistent.
"""
function _target_inputs(inputs::NamedTuple, TA, DIC, temp_c, pres_bar)
    cleared = NamedTuple{DERIVED_PARAMETERS}(ntuple(_ -> nothing, length(DERIVED_PARAMETERS)))

    return merge(inputs, cleared, (
        temp_c = something(temp_c, inputs.temp_c),
        pres_bar = something(pres_bar, inputs.pres_bar),
        TA = TA,
        DIC = DIC,
    ))
end

"`result`'s inputs, re-targeted, with TA and DIC taken from its solved state."
_target_inputs(result::CarbonateResult, temp_c, pres_bar) =
    _target_inputs(result.inputs, result.val.TA, result.val.DIC, temp_c, pres_bar)

"""
Whether a result was solved from TA and DIC themselves.

Then both are already known at any conditions, and re-solving needs nothing from the original
solve: the new state follows from TA, DIC and the new conditions alone.
"""
_solved_from_TA_DIC(result::CarbonateResult) = issetequal(result.input_keys, (:TA, :DIC))

"""
    adjust(result; temp_c, pres_bar, σ_temp_c, σ_pres_bar)
    adjust(result, temp_c, pres_bar = nothing, σ_temp_c = nothing, σ_pres_bar = nothing)

Re-solve a carbonate system at a different temperature and/or pressure.

`result` comes from [`carbon_system`](@ref), [`whole_system`](@ref) or a
[`CarbonateSystem`](@ref) whose scope includes `:carbon`. This returns the state of the same
water — the same TA, DIC, totals and composition — at the new conditions. Arguments left as
`nothing` keep their original value.

A common use is reporting samples measured in the laboratory at the temperature and pressure
they were collected at; others include normalising a dataset to one temperature, or following
a water mass as it warms or sinks. The name follows PyCO2SYS v2, whose equivalent is also
`adjust`.

Salinity is not an argument: water with a different salinity is different water, not the same
water at other conditions.

# Examples
```jldoctest
julia> measured = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0);   # on deck

julia> in_situ = adjust(measured, temp_c = 2.0, pres_bar = 400.0);

julia> in_situ.pHtot
8.254156644905194
```

# Arguments

- `result`: a solved state from [`carbon_system`](@ref), [`whole_system`](@ref) or a
  `CarbonateSystem`.
- `temp_c = nothing`: the new temperature, °C. Keeps the original value when `nothing`.
- `pres_bar = nothing`: the new pressure, bar (hydrostatic; 0 at the surface).
- `σ_temp_c`, `σ_pres_bar = nothing`: uncertainty in the new conditions, which are often less
  well known than the original ones (an in-situ temperature, say).

# Extended help

Everything needed to re-solve — the conservative totals, the seawater composition, and every
setting — is carried on `result`, so nothing has to be restated.

**What gets recalculated.** Without uncertainties this is one solve, at the new conditions,
from the TA and DIC already on `result`. With uncertainties, how much is recalculated depends
on what the original calculation was given:

- **From TA and DIC:** the same single solve. Both are conserved, so the original conditions
  cannot affect the answer, and there is nothing earlier to differentiate through.
- **From anything else** (pH and TA, say): the original solve is repeated inside the
  differentiation, because TA or DIC was *derived* there using constants at the original
  conditions, and its uncertainty — including any from those conditions — has to reach the
  answer.

If you have TA and DIC and only want the system at one set of conditions, you need no
`adjust` at all: solve directly at those conditions, `carbon_system(TA = …, DIC = …,
temp_c = …, pres_bar = …)`.

`result` must come straight from a calculation, not from `adjust` itself. To report the same
water at several sets of conditions, adjust the original result once for each.

The second form takes the same arguments positionally, which is what makes it broadcastable —
Julia will not broadcast over keyword arguments. Scalars broadcast against vectors, so a value
shared by every sample is just a number.

Uncertainties given to the original calculation are carried on `result` and do not have to be
restated. A calculation carrying uncertainties is built with `varying_errors` and called
positionally — the presets take no `errors` argument:

```jldoctest adjust
julia> measured = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC),
                                  TA = 2300.0, DIC = 2000.0, temp_c = 20.0)(2.0, 2.0);

julia> adjust(measured; temp_c = 2.0, σ_temp_c = 0.5).err.pHtot
0.009597793412851123
```

The two contributions are independent and combine in quadrature. Whether uncertainty in the
*original* temperature reaches the new state at all depends on what was measured: with TA and
DIC it does not, because both are conservative and carried across unchanged, but with pH and
DIC it does, because TA is derived using constants at the original temperature.

A whole cast at once, carrying per-sample uncertainty in the in-situ temperature:

```jldoctest adjust
julia> in_situ = adjust.([measured, measured], [2.0, 4.0], 400.0, [0.5, 0.2]);

julia> getproperty.(in_situ, :pHtot)
2-element Vector{Float64}:
 8.254156644905194
 8.223138522847922
```

See also [`carbon_system`](@ref), [`CarbonateResult`](@ref).
"""
function adjust(result::CarbonateResult;
                temp_c = nothing,
                pres_bar = nothing,
                σ_temp_c = nothing,
                σ_pres_bar = nothing,
                unsupported...)
    isempty(unsupported) || _reject_adjust_keywords(keys(unsupported))
    return _adjust(result, temp_c, pres_bar, σ_temp_c, σ_pres_bar)
end

"""
The positional form, which is the one that broadcasts.

`temp_c` has no default, unlike every argument after it. If it had one this method would also
define the single-argument call and overwrite the keyword method above, which already covers it.
Pass `nothing` explicitly to hold a condition at its original value while varying a later one:
`adjust.(results, nothing, nothing, σ)`.
"""
adjust(result::CarbonateResult, temp_c, pres_bar = nothing,
       σ_temp_c = nothing, σ_pres_bar = nothing) =
    _adjust(result, temp_c, pres_bar, σ_temp_c, σ_pres_bar)

"The one path both call forms take."
function _adjust(result::CarbonateResult, temp_c, pres_bar, σ_temp_c, σ_pres_bar)
    _reject_vector_arguments(temp_c, pres_bar, σ_temp_c, σ_pres_bar)

    # Adjusting an adjusted result would chain uncertainty through two condition changes when
    # the water only ever had one original state. `source` is set only by this function, so
    # having one is what marks a result as already adjusted.
    isnothing(result.source) || throw(ArgumentError(
        "this result was produced by `adjust`, which only takes a result straight from a " *
        "calculation.\nTo report the same water at other conditions, adjust the original " *
        "result again rather than this one."))

    :carbon in result.system || throw(ArgumentError(
        "adjust carries TA and DIC to the new conditions, so it needs a result whose scope " *
        "includes :carbon; this one was solved with scope $(result.system)."))

    core = _core_for(result.system)
    target = _target_inputs(result, temp_c, pres_bar)

    original = result.inputs
    resolved_temp_c = something(temp_c, original.temp_c)
    resolved_pres_bar = something(pres_bar, original.pres_bar)

    # Both exits below call the core directly rather than going through `_run`, so this is the
    # only place the new conditions are checked at all. Without it a target of 200 °C or a
    # negative pressure returns a number in silence, where a direct calculation would warn or
    # throw.
    #
    # Here rather than inside `composed`: `_check_conditions` returns early for anything that is
    # not `Real`, so on the AD path it would meet `Dual`s and never fire. At this point the
    # values are still plain numbers.
    #
    # `sal` comes from the original because it is not a degree of freedom here, but it still
    # reaches the equilibrium constants at the new conditions, so it still needs checking.
    _check_conditions(resolved_temp_c, original.sal, resolved_pres_bar)

    # No `_check_determinacy`: `_target_inputs` clears `DERIVED_PARAMETERS` and sets TA and DIC
    # from the solved state, so the target is determinate by construction rather than by luck.

    condition_errors = _target_errors(σ_temp_c, σ_pres_bar)
    ad_errors = merge(something(result.input_errors, NamedTuple()), condition_errors)

    # What this state was adjusted from. `inputs` is about to be overwritten with the new
    # conditions and the derived parameters cleared, and `condition_errors` is recorded nowhere
    # else, so without this an adjusted state cannot reconstruct its own chain — which is what
    # propagating an uncertainty through anything computed *from* it needs.
    source = (inputs = original, errors = result.input_errors,
              condition_errors = condition_errors)

    if isempty(ad_errors)
        return CarbonateResult(core(; target...), nothing, result.input_keys,
                               _settings(target), target, result.system, nothing;
                               source = source)
    end

    # Differentiate the whole chain, from the original independent measurements through to
    # the new state, rather than propagating the original output uncertainties forward.
    #
    # That distinction matters. Salinity, for instance, affects the equilibrium constants at
    # *both* conditions, so carrying σ(TA)/σ(DIC) forward and re-propagating σ(S) here would
    # count salinity twice whenever TA or DIC was itself derived - and would lose the
    # correlation between them. Differentiating from the independent inputs gets every such
    # case right without having to enumerate them. Viable only because the constants are
    # pure Julia and AD-transparent.
    #
    # The new conditions enter under distinct names so they cannot collide with the original
    # ones.
    ad_inputs = merge(original, (target_temp_c = resolved_temp_c,
                                 target_pres_bar = resolved_pres_bar))
    from_TA_DIC = _solved_from_TA_DIC(result)

    function composed(; kwargs...)
        kw = NamedTuple(kwargs)
        original_kw = Base.structdiff(kw, NamedTuple{(:target_temp_c, :target_pres_bar)})

        # Measured TA and DIC are already the conservative pair, so the original solve would
        # only hand them back: go straight to the new conditions. The derivatives are the same
        # either way, since nothing at the original conditions can reach TA or DIC.
        from_TA_DIC && return core(; _target_inputs(original_kw, kw.TA, kw.DIC,
                                                    kw.target_temp_c, kw.target_pres_bar)...)

        # Otherwise TA or DIC was derived at the original conditions, so that solve is part
        # of the chain and has to be differentiated through.
        solved = core(; original_kw...)
        return core(; _target_inputs(original_kw, solved.TA, solved.DIC,
                                     kw.target_temp_c, kw.target_pres_bar)...)
    end

    res_data = propagate_errors(composed; inputs=ad_inputs, errors=ad_errors)
    return CarbonateResult(res_data.val, res_data.err, result.input_keys,
                           _settings(target), target, result.system, result.input_errors;
                           source = source)
end

"""
    at_collection_conditions(result; temp_c, pres_bar, σ_temp_c, σ_pres_bar)

Deprecated: the former name of [`adjust`](@ref), which takes the same arguments.

Kept so existing scripts keep working. It forwards to `adjust`, and warns when Julia is run
with deprecation warnings on (as `Pkg.test` does).
"""
function at_collection_conditions(result::CarbonateResult, args...; kwargs...)
    Base.depwarn("`at_collection_conditions` is deprecated; it has been renamed `adjust`, " *
                 "which takes the same arguments.", :at_collection_conditions)
    return adjust(result, args...; kwargs...)
end

"""
The uncertainties in the new conditions, under the internal names the AD path uses.

One method per combination rather than a branch, so which uncertainties are present is settled
by the argument types and folds away on specialisation — a broadcast passes `Float64` where a σ
column was given and `Nothing` where a slot was skipped, and each call site compiles to one
concrete NamedTuple. Building the key set at runtime instead costs an allocation and a dynamic
dispatch on every call (see `CarbonateSystem`).

The `target_` prefix is applied here rather than mapped over user-supplied names, which is what
makes a collision with an original uncertainty of the same name impossible.
"""
_target_errors(::Nothing, ::Nothing) = NamedTuple()
_target_errors(σ_temp_c, ::Nothing) = (target_temp_c = σ_temp_c,)
_target_errors(::Nothing, σ_pres_bar) = (target_pres_bar = σ_pres_bar,)
_target_errors(σ_temp_c, σ_pres_bar) = (target_temp_c = σ_temp_c,
                                        target_pres_bar = σ_pres_bar)

"""
Reject a vector handed to an argument that describes one sample.

`adjust.(results; temp_c = temps)` is syntactically valid and broadcasts over `results` alone,
handing the whole vector of temperatures to every call. Unchecked it surfaces as
`isless(::Int64, ::Vector)` from inside the solver, or `convert(Float64, ::Vector)` from the
uncertainty path — neither of which names the mistake.

Written out rather than looped over, so each check folds away for the concrete argument types a
broadcast produces.
"""
function _reject_vector_arguments(temp_c, pres_bar, σ_temp_c, σ_pres_bar)
    temp_c     isa AbstractArray && _vector_argument_error(:temp_c, temp_c)
    pres_bar   isa AbstractArray && _vector_argument_error(:pres_bar, pres_bar)
    σ_temp_c   isa AbstractArray && _vector_argument_error(:σ_temp_c, σ_temp_c)
    σ_pres_bar isa AbstractArray && _vector_argument_error(:σ_pres_bar, σ_pres_bar)
    return nothing
end

@noinline _vector_argument_error(name, value) = throw(ArgumentError(
    "$name was given a $(typeof(value)), but one call solves one sample.\n" *
    "Broadcast the positional form instead: adjust.(results, temperatures)"))

"""
Reject keyword arguments `adjust` does not take.

Named individually for the two people will actually reach for: `sal`, because a different
salinity is a reasonable thing to want and a wrong thing to ask for here, and `errors`, because
uncertainty in a measured quantity belongs to the original calculation and is carried here on
the result.
"""
function _reject_adjust_keywords(names)
    lines = map(collect(names)) do name
        if name === :sal || name === :σ_sal
            "$name: water with a different salinity is different water, not the same water at " *
            "other conditions, so adjust does not change it.\n    If the salinity genuinely " *
            "differs, that is a direct solve — broadcast a CarbonateSystem with :sal in " *
            "`varying`.\n    For its uncertainty, rebuild the original calculation with " *
            "CarbonateSystem(:carbon; varying_errors = (:sal,))."
        elseif name === :errors
            "errors: give uncertainties in the new conditions as σ_temp_c and σ_pres_bar.\n" *
            "    Uncertainty in a measured quantity belongs to the original calculation, so it " *
            "goes in the solver's `varying_errors` and is carried here on the result " *
            "automatically."
        else
            "$name is not something adjust changes; expected temp_c, pres_bar, σ_temp_c " *
            "or σ_pres_bar"
        end
    end
    throw(ArgumentError("unrecognised argument(s):\n  " * join(lines, "\n  ")))
end
