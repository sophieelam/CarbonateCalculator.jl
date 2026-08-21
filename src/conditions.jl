# Stage 2 of the two-stage conditions API.
#
# Stage 1 (`carbon_system` / `whole_system`) solves the system at the conditions the sample
# was *measured* at. This file re-solves that result at the conditions the sample was
# *collected* at.
#
# Salinity is not a parameter: a water sample's salinity does not change between collection
# and measurement, only its temperature and pressure do.

"""
Carbonate parameters that describe the *solved state* rather than the sample itself. They
must be cleared before re-solving, or the target calculation would be driven by the values
measured at the original conditions.
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
Build the argument list that re-solves `result` at the target conditions.

TA and DIC are taken from the *solved* state rather than the original inputs, because they
are the conservative pair and are what carries the chemistry across the condition change.
Everything else — the totals, the seawater composition, and every `*_method` — is inherited
from `result.inputs`, which is what keeps the two conditions consistent with each other.
"""
function _target_inputs(result::CarbonateResult, temp_c, pres_bar)
    inputs = result.inputs
    cleared = NamedTuple{DERIVED_PARAMETERS}(ntuple(_ -> nothing, length(DERIVED_PARAMETERS)))

    return merge(inputs, cleared, (
        temp_c = something(temp_c, inputs.temp_c),
        pres_bar = something(pres_bar, inputs.pres_bar),
        TA = result.val.TA,
        DIC = result.val.DIC,
    ))
end

"""
    at_collection_conditions(result; temp_c, pres_bar, σ_temp_c, σ_pres_bar)
    at_collection_conditions(result, temp_c, pres_bar = nothing,
                             σ_temp_c = nothing, σ_pres_bar = nothing)

Re-solve a carbonate system at a different temperature and/or pressure.

`result` comes from [`carbon_system`](@ref) or [`whole_system`](@ref) and describes the
system at the conditions it was *measured* at. This returns the state of the same water at
the conditions it was *collected* at. Arguments left as `nothing` keep their measured value.

Salinity is not an argument: a sample's salinity does not change between collection and
measurement.

# Examples
```jldoctest
julia> measured = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0);   # on deck

julia> collected = at_collection_conditions(measured, temp_c = 2.0, pres_bar = 400.0);

julia> collected.pHtot
8.254156644905194
```

# Arguments

- `result`: a measured state from [`carbon_system`](@ref) or [`whole_system`](@ref).
- `temp_c = nothing`: collection temperature, °C. Keeps the measured value when `nothing`.
- `pres_bar = nothing`: collection pressure, bar.
- `σ_temp_c`, `σ_pres_bar = nothing`: uncertainty in the collection conditions, which are
  often less well known than the conditions the sample was measured at.

# Extended help

Everything needed to re-solve — the conservative totals, the seawater composition, and every
method choice — is carried on `result`, so nothing has to be restated.

`result` must be a measured state. A sample has one set of collection conditions, so the
result of this function cannot be passed back into it; to report the same measurement at some
other conditions, carry the measured result again.

The second form takes the same arguments positionally, which is what makes it broadcastable —
Julia will not broadcast over keyword arguments. Scalars broadcast against vectors, so a value
shared by every sample is just a number.

Uncertainties given to the measurement are carried on `result` and do not have to be restated.
A measurement carrying uncertainties is built with `varying_errors` and called positionally —
the presets take no `errors` argument:

```jldoctest collection
julia> measured = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC),
                                  TA = 2300.0, DIC = 2000.0, temp_c = 20.0)(2.0, 2.0);

julia> at_collection_conditions(measured; temp_c = 2.0, σ_temp_c = 0.5).err.pHtot
0.009597793412851123
```

The two contributions are independent and combine in quadrature. Whether uncertainty in the
*measurement* temperature reaches the collection state at all depends on what was measured:
with TA and DIC it does not, because both are conservative and carried across unchanged, but
with pH and DIC it does, because TA is derived using constants at the measurement temperature.

A whole cast at once, carrying per-sample uncertainty in the collection temperature:

```jldoctest collection
julia> collected = at_collection_conditions.([measured, measured], [2.0, 4.0],
                                             400.0, [0.5, 0.2]);

julia> getproperty.(collected, :pHtot)
2-element Vector{Float64}:
 8.254156644905194
 8.223138522847922
```

See also [`carbon_system`](@ref), [`CarbonateResult`](@ref).
"""
function at_collection_conditions(result::CarbonateResult;
                                  temp_c = nothing,
                                  pres_bar = nothing,
                                  σ_temp_c = nothing,
                                  σ_pres_bar = nothing,
                                  unsupported...)
    isempty(unsupported) || _reject_collection_keywords(keys(unsupported))
    return _at_collection_conditions(result, temp_c, pres_bar, σ_temp_c, σ_pres_bar)
end

"""
The positional form, which is the one that broadcasts.

`temp_c` has no default, unlike every argument after it. If it had one this method would also
define the single-argument call and overwrite the keyword method above, which already covers it.
Pass `nothing` explicitly to hold a condition at its measured value while varying a later one:
`at_collection_conditions.(results, nothing, nothing, σ)`.
"""
at_collection_conditions(result::CarbonateResult, temp_c, pres_bar = nothing,
                         σ_temp_c = nothing, σ_pres_bar = nothing) =
    _at_collection_conditions(result, temp_c, pres_bar, σ_temp_c, σ_pres_bar)

"The one path both call forms take."
function _at_collection_conditions(result::CarbonateResult, temp_c, pres_bar,
                                   σ_temp_c, σ_pres_bar)
    _reject_vector_arguments(temp_c, pres_bar, σ_temp_c, σ_pres_bar)

    # A sample has one set of collection conditions, so a second carry has nothing to mean.
    # Carrying a result is the only thing that gives it a `source`, so having one is what marks
    # it as already carried.
    isnothing(result.source) || throw(ArgumentError(
        "this result is already at its collection conditions, and a sample has only one set " *
        "of those.\nTo report the same measurement somewhere else, carry the measured result " *
        "again rather than this one."))

    core = _core_for(result.system)
    target = _target_inputs(result, temp_c, pres_bar)

    measured = result.inputs
    resolved_temp_c = something(temp_c, measured.temp_c)
    resolved_pres_bar = something(pres_bar, measured.pres_bar)

    # Both exits below call the core directly rather than going through `_run`, so this is the
    # only place the target conditions are checked at all. Without it a target of 200 °C or a
    # negative pressure returns a number in silence, where stage 1 would warn or throw.
    #
    # Here rather than inside `composed`: `_check_conditions` returns early for anything that is
    # not `Real`, so on the AD path it would meet `Dual`s and never fire. At this point the
    # values are still plain numbers.
    #
    # `sal` comes from the measurement because it is not a degree of freedom here, but it still
    # reaches the equilibrium constants at the target conditions, so it still needs checking.
    _check_conditions(resolved_temp_c, measured.sal, resolved_pres_bar)

    # No `_check_determinacy`: `_target_inputs` clears `DERIVED_PARAMETERS` and sets TA and DIC
    # from the solved state, so the target is determinate by construction rather than by luck.

    condition_errors = _target_errors(σ_temp_c, σ_pres_bar)
    ad_errors = merge(something(result.input_errors, NamedTuple()), condition_errors)

    # What this state was carried from. `inputs` is about to be overwritten with the collection
    # conditions and the derived parameters cleared, and `condition_errors` is recorded nowhere
    # else, so without this a collection state cannot reconstruct its own chain — which is what
    # propagating an uncertainty through anything computed *from* it needs.
    source = (inputs = measured, errors = result.input_errors,
              condition_errors = condition_errors)

    if isempty(ad_errors)
        return CarbonateResult(core(; target...), nothing, result.input_keys,
                               _settings(target), target, result.system, nothing;
                               source = source)
    end

    # Differentiate the whole chain, from the original independent measurements through to
    # the target state, rather than propagating stage 1's output uncertainties forward.
    #
    # That distinction matters. Salinity, for instance, affects the equilibrium constants at
    # *both* conditions, so carrying σ(TA)/σ(DIC) forward and re-propagating σ(S) here would
    # count salinity twice whenever TA or DIC was itself derived - and would lose the
    # correlation between them. Differentiating from the independent inputs gets every such
    # case right without having to enumerate them. Viable only because the constants are
    # pure Julia and AD-transparent.
    #
    # The target conditions enter under distinct names so they cannot collide with the
    # measured ones.
    ad_inputs = merge(measured, (target_temp_c = resolved_temp_c,
                                 target_pres_bar = resolved_pres_bar))

    function composed(; kwargs...)
        kw = NamedTuple(kwargs)
        stage1_kw = Base.structdiff(kw, NamedTuple{(:target_temp_c, :target_pres_bar)})
        stage1 = core(; stage1_kw...)
        measured_result = CarbonateResult(stage1, nothing, result.input_keys,
                                          _settings(stage1_kw), stage1_kw, result.system,
                                          nothing)
        return core(; _target_inputs(measured_result, kw.target_temp_c, kw.target_pres_bar)...)
    end

    res_data = propagate_errors(composed; inputs=ad_inputs, errors=ad_errors)
    return CarbonateResult(res_data.val, res_data.err, result.input_keys,
                           _settings(target), target, result.system, result.input_errors;
                           source = source)
end

"""
The collection-condition uncertainties, under the internal names the AD path uses.

One method per combination rather than a branch, so which uncertainties are present is settled
by the argument types and folds away on specialisation — a broadcast passes `Float64` where a σ
column was given and `Nothing` where a slot was skipped, and each call site compiles to one
concrete NamedTuple. Building the key set at runtime instead costs an allocation and a dynamic
dispatch on every call (see `CarbonateSystem`).

The `target_` prefix is applied here rather than mapped over user-supplied names, which is what
makes a collision with a stage-1 uncertainty of the same name impossible.
"""
_target_errors(::Nothing, ::Nothing) = NamedTuple()
_target_errors(σ_temp_c, ::Nothing) = (target_temp_c = σ_temp_c,)
_target_errors(::Nothing, σ_pres_bar) = (target_pres_bar = σ_pres_bar,)
_target_errors(σ_temp_c, σ_pres_bar) = (target_temp_c = σ_temp_c,
                                        target_pres_bar = σ_pres_bar)

"""
Reject a vector handed to an argument that describes one sample.

`at_collection_conditions.(results; temp_c = temps)` is syntactically valid and broadcasts over
`results` alone, handing the whole vector of temperatures to every call. Unchecked it surfaces
as `isless(::Int64, ::Vector)` from inside the solver, or `convert(Float64, ::Vector)` from the
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
    "Broadcast the positional form instead: at_collection_conditions.(results, temperatures)"))

"""
Reject keyword arguments that are not collection conditions.

Named individually for the two people will actually reach for: `sal`, because a per-sample
collection salinity is a reasonable thing to want and a wrong thing to ask for here, and
`errors`, because uncertainty in a measured quantity belongs to the measurement and is carried
here on the result.
"""
function _reject_collection_keywords(names)
    lines = map(collect(names)) do name
        if name === :sal || name === :σ_sal
            "$name: a sample's salinity does not change between collection and measurement, " *
            "so it is not a collection condition.\n    If the collection salinity genuinely " *
            "differs then it is not the same water, and that is a direct solve — broadcast a " *
            "CarbonateSystem with :sal in `varying`.\n    For its uncertainty, rebuild the " *
            "measurement with CarbonateSystem(:carbon; varying_errors = (:sal,))."
        elseif name === :errors
            "errors: give collection-condition uncertainties as σ_temp_c and σ_pres_bar.\n" *
            "    Uncertainty in a measured quantity belongs to the measurement, so it goes in " *
            "the solver's `varying_errors` and is carried here on the result automatically."
        else
            "$name is not a collection condition; expected temp_c, pres_bar, σ_temp_c " *
            "or σ_pres_bar"
        end
    end
    throw(ArgumentError("unrecognised argument(s):\n  " * join(lines, "\n  ")))
end
