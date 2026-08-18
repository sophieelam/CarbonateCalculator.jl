# Stage 2 of the two-stage conditions API.
#
# Stage 1 (`carbon_system` / `whole_system`) solves the system at the conditions the sample
# was *measured* at. This file re-solves that result at the conditions the sample was
# *collected* at.
#
# Salinity is deliberately not a parameter: a water sample's salinity does not change
# between collection and measurement, only its temperature and pressure do.

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

A result records the *scope* it was solved with, so re-solving it needs no guesswork about
which entry point it came from — and works for any scope, not just the two that used to have
their own cores.
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
    recalculate_at_target_conditions(result; temp_c=nothing, pres_bar=nothing, errors=nothing)

Re-solve a carbonate system at a different temperature and/or pressure.

`result` comes from [`carbon_system`](@ref) or [`whole_system`](@ref) and describes the
system at the conditions it was *measured* at. This returns the state of the same water at
the conditions it was *collected* at. Arguments left as `nothing` keep their measured value.

Salinity is not an argument: a sample's salinity does not change between collection and
measurement.

Everything needed to re-solve — the conservative totals, the seawater composition, and every
method choice — is carried on `result`, so nothing has to be restated.

# Examples
```julia
measured = carbon_system(TA=2300.0, DIC=2000.0, temp_c=20.0)   # on deck
collected = recalculate_at_target_conditions(measured, temp_c=2.0, pres_bar=400.0)
collected.pHtot
```

Uncertainties given to stage 1 are carried on `result` and do not have to be restated.
`errors` here adds uncertainty in the *target conditions*, which are often less well known
than the conditions the sample was measured at:

```julia
measured = carbon_system(TA=2300.0, DIC=2000.0, temp_c=20.0, errors=(TA=2.0, DIC=2.0))
recalculate_at_target_conditions(measured; temp_c=2.0, errors=(temp_c=0.5,)).err.pHtot
```
"""
function recalculate_at_target_conditions(result::CarbonateResult;
                                          temp_c=nothing,
                                          pres_bar=nothing,
                                          errors=nothing)
    core = _core_for(result.system)
    target = _target_inputs(result, temp_c, pres_bar)

    measured = result.inputs
    resolved_temp_c = something(temp_c, measured.temp_c)
    resolved_pres_bar = something(pres_bar, measured.pres_bar)

    ad_errors = merge(something(result.input_errors, NamedTuple()),
                      isnothing(errors) ? NamedTuple() : _rename_target_errors(errors))

    if isempty(ad_errors)
        return CarbonateResult(core(; target...), nothing, result.input_keys,
                               _settings(target), target, result.system, nothing)
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
                           _settings(target), target, result.system, result.input_errors)
end

"""
Map user-facing `temp_c`/`pres_bar` uncertainties onto the internal target-condition names,
so they cannot be confused with uncertainty in the conditions the sample was measured at.
"""
function _rename_target_errors(errors::NamedTuple)
    renamed = map(keys(errors)) do k
        k === :temp_c   ? :target_temp_c   :
        k === :pres_bar ? :target_pres_bar : k
    end
    return NamedTuple{Tuple(renamed)}(values(errors))
end
