# Buffer factors: how one carbonate variable responds to another.
#
# The Revelle factor is one instance of a general quantity — ∂X/∂Y at constant Z — so rather
# than hardcoding each, the system is re-solved with (Y, Z) as its defining pair and
# differentiated. The constants are AD-transparent, so this is exact rather than a finite
# difference.
#
# These are deliberately not part of a solve. Each one costs a re-solve in Dual arithmetic,
# most callers want none of them, and there are more of them than anyone would want by default.

"""
Every parameter that constrains the carbonate system, from the groups the determinacy check
uses. All of them are stripped before the chosen pair goes back in, or whatever the result was
originally solved from would over-determine it.

Taken from `PARAMETER_GROUPS` rather than listed again here, so a new parameter cannot be
added to the solver and forgotten here — and so the boron and isotope constraints are covered
without naming them.
"""
const CONSTRAINING_PARAMETERS = Tuple(unique(Iterators.flatten(values(PARAMETER_GROUPS))))

# Built once. Rebuilding a NamedTuple from a key set per call is the expensive pattern in this
# package, and this one never varies.
const _CLEARED_CONSTRAINTS =
    NamedTuple{CONSTRAINING_PARAMETERS}(ntuple(_ -> nothing, length(CONSTRAINING_PARAMETERS)))

"""
The inputs that re-solve `result` with every constraint stripped but `constant`, ready for the
denominator to be substituted in.

`constant` is taken from the *solved* state rather than the original inputs, so it can be a
quantity the result derived rather than one it was given — holding pH constant is legitimate
whether or not pH was measured.
"""
function _gradient_inputs(result::CarbonateResult, constant::Symbol)
    held = NamedTuple{(constant,)}((getproperty(result.val, constant),))
    return merge(result.inputs, _CLEARED_CONSTRAINTS, held)
end

"""
    calc_gradient(result, numerator, denominator; constant)

∂`numerator`/∂`denominator`, holding `constant` fixed, by automatic differentiation.

`(denominator, constant)` becomes the pair that defines the system while the derivative is
taken, so it has to determine it — the same rule the solver applies to its own inputs, and
checked the same way. Everything else about the sample is inherited from `result`.

What is held fixed is half the definition: ∂pCO₂/∂DIC at constant TA is the Revelle factor,
while ∂pCO₂/∂DIC at constant pH is a different quantity roughly ten times smaller.

Temperature, salinity and pressure are held by construction — they are not carbonate
parameters, so nothing perturbs them.

```julia
result = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0)
calc_gradient(result, :pCO₂, :DIC; constant = :TA)     # ∂pCO₂/∂DIC
calc_gradient(result, :TA, :pHtot; constant = :DIC)    # buffer capacity, ∂TA/∂pH
```

See [`calc_relative_gradient`](@ref) for the normalised form, and [`revelle_factor`](@ref).
"""
function calc_gradient(result::CarbonateResult, numerator::Symbol, denominator::Symbol;
                       constant::Symbol)
    base, start = _validated_request(result, denominator, constant)
    return _derivative_at(result.system, base, numerator, denominator, start)
end

"""
∂`numerator`/∂`denominator` for a system already pinned by `base`, evaluated at `start`.

`base` carries the held constant and everything about the sample that is not a carbonate
constraint, so substituting a denominator into it names one system.
"""
function _derivative_at(scope, base, numerator::Symbol, denominator::Symbol, start)
    solved_at(y) = getproperty(_solve_core(scope, merge(base, NamedTuple{(denominator,)}((y,)))),
                               numerator)
    return ForwardDiff.derivative(solved_at, start)
end

"""
Check a gradient request and return what computing it needs.

The determinacy check is the solver's own, so an ill-posed pair is refused for the same reason
and with the same message it would get as input — including that pCO₂ and fCO₂ are one quantity
rather than two independent ones.
"""
function _validated_request(result::CarbonateResult, denominator::Symbol, constant::Symbol)
    denominator === constant && throw(ArgumentError(
        "$denominator cannot be both varied and held constant"))

    base = _gradient_inputs(result, constant)
    start = getproperty(result.val, denominator)

    _check_determinacy(merge(base, NamedTuple{(denominator,)}((start,))), result.system;
                       require_two = result.system === (:carbon,))

    return base, start
end

"""
    calc_relative_gradient(result, numerator, denominator; constant)

The normalised gradient, `(∂X/∂Y)·(Y/X)` — a fractional change in one quantity per fractional
change in another, and therefore dimensionless.

This is the form the named buffer factors take, the Revelle factor among them. Arguments are
otherwise exactly [`calc_gradient`](@ref)'s.

```julia
result = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0)
calc_relative_gradient(result, :pCO₂, :DIC; constant = :TA)   # the Revelle factor
```
"""
function calc_relative_gradient(result::CarbonateResult, numerator::Symbol,
                                denominator::Symbol; constant::Symbol)
    gradient = calc_gradient(result, numerator, denominator; constant = constant)
    return gradient * getproperty(result.val, denominator) / getproperty(result.val, numerator)
end

"""
    calc_gradient_uncertainty(result, numerator, denominator; constant, relative = false)

The uncertainty in a gradient, propagated from the measurement uncertainties `result` carries.

Returns `nothing` when the result has none, since there is then nothing to propagate.

A gradient is itself a derivative, so its uncertainty is a second derivative — the whole chain
is differentiated twice over, once inside to take the gradient and once outside to see how it
moves with each measured input.

**The evaluation point moves with the inputs, and it has to.** A gradient is taken at the
sample's own state, and both the point it is evaluated at and the value held constant are
functions of what was measured. Freezing them at their solved values and differentiating only
the surrounding solve returns **exactly zero** for uncertainty in any carbonate parameter,
because `_gradient_inputs` strips those before re-solving and the perturbation never reaches
the calculation. Worse, it still returns a plausible number for σ on temperature or salinity,
which survive the stripping — so the mistake is invisible in the easiest thing to test with.
Re-solving from the measured inputs is the same reasoning that makes `at_collection_conditions`
differentiate its whole chain rather than propagate an intermediate forward.

For a result carried to its collection conditions the chain runs through both stages, from the
measured inputs to the collection state to the gradient, using the provenance
`at_collection_conditions` records on `result.source`. Uncertainty in the collection conditions
themselves is included alongside the measurement's. Taking the uncertainties at face value at the
collection conditions instead would be wrong whenever the measured pair was not the conservative
one: `input_errors` may name `pHtot`, which the collection state does not carry.

Uses the same first-order, uncorrelated model as [`propagate_errors`](@ref), and inherits its
one assumption.

```julia
measured = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC),
                           TA = 2300.0, DIC = 2000.0, temp_c = 20.0)(2.0, 2.0)
calc_gradient_uncertainty(measured, :pCO₂, :DIC; constant = :TA)
```
"""
function calc_gradient_uncertainty(result::CarbonateResult, numerator::Symbol,
                                   denominator::Symbol; constant::Symbol,
                                   relative::Bool = false)
    _validated_request(result, denominator, constant)

    inputs, errors = _uncertainty_chain(result)
    (isnothing(errors) || isempty(errors)) && return nothing

    scope = result.system
    names = Val(keys(errors))
    solve_to_state = _state_solver(result)

    function gradient_from(perturbations)
        perturbed = _perturbed(inputs, names, perturbations)

        # Re-solved, so the point the gradient is taken at and the value held constant both
        # carry the perturbation. This is what freezing them would throw away.
        state, base_inputs = solve_to_state(perturbed)
        start = getproperty(state, denominator)
        base = merge(base_inputs, _CLEARED_CONSTRAINTS,
                     NamedTuple{(constant,)}((getproperty(state, constant),)))

        gradient = _derivative_at(scope, base, numerator, denominator, start)
        relative || return gradient

        # Normalised here rather than outside, because the numerator and denominator it is
        # scaled by are uncertain too — scaling σ afterwards would miss that.
        return gradient * start / getproperty(state, numerator)
    end

    measured = Float64[getproperty(inputs, name) for name in keys(errors)]
    sensitivities = ForwardDiff.gradient(gradient_from, measured)
    sigmas = Float64[getproperty(errors, name) for name in keys(errors)]

    return sqrt(sum((sensitivities .* sigmas) .^ 2))
end

"""
The inputs a gradient's uncertainty is differentiated with respect to, and their σ.

A measurement is differentiated with respect to what it was given. A collection state is
differentiated with respect to what the *measurement* was given, plus the collection conditions
themselves — its own `inputs` are the target ones, where a pH the measurement supplied has been
cleared and the measurement's temperature overwritten, so perturbing those would be perturbing
the wrong thing.
"""
function _uncertainty_chain(result::CarbonateResult)
    source = result.source
    isnothing(source) && return result.inputs, result.input_errors

    conditions = (target_temp_c = result.inputs.temp_c,
                  target_pres_bar = result.inputs.pres_bar)
    return merge(source.inputs, conditions),
           merge(something(source.errors, NamedTuple()), source.condition_errors)
end

"""
How to get from a perturbed input set to the state the gradient is taken at, and to the inputs
that name that state.

For a measurement that is one solve. For a collection state it is the same two-stage chain
`at_collection_conditions` differentiates: solve at the measurement, carry the conservative pair
across to the collection conditions, solve again.
"""
function _state_solver(result::CarbonateResult)
    scope = result.system

    isnothing(result.source) &&
        return perturbed -> (_solve_core(scope, perturbed), perturbed)

    cleared = NamedTuple{DERIVED_PARAMETERS}(ntuple(_ -> nothing, length(DERIVED_PARAMETERS)))

    return function (perturbed)
        measured_inputs = Base.structdiff(perturbed,
                                          NamedTuple{(:target_temp_c, :target_pres_bar)})
        measured_state = _solve_core(scope, measured_inputs)

        target = merge(measured_inputs, cleared,
                       (temp_c = perturbed.target_temp_c,
                        pres_bar = perturbed.target_pres_bar,
                        TA = measured_state.TA, DIC = measured_state.DIC))

        return _solve_core(scope, target), target
    end
end

"""
    revelle_factor(result)

The Revelle factor: the fractional change in CO₂ per fractional change in DIC, at constant
alkalinity.

A preset over [`calc_relative_gradient`](@ref). `fCO₂` and `pCO₂` differ by a factor fixed at
constant temperature, salinity and pressure, so either gives the same number here.

```julia
revelle_factor(carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0))
```
"""
revelle_factor(result::CarbonateResult) =
    calc_relative_gradient(result, :fCO₂, :DIC; constant = :TA)

"""
    with_gradient(result, name, numerator, denominator; constant, relative = false)

`result` with a gradient added to its computed values under `name`, so it travels with the
sample and lands in a column of its own when a vector of results goes into a table.

The name is given rather than derived from the arguments: a generated `∂pCO₂_∂DIC` would have
to be built with a runtime `Symbol`, which is the expensive pattern this package avoids, and
would read worse than whatever the caller would have called it.

A new result is returned — `CarbonateResult` is immutable, and its values are concretely typed,
which is what the Tables.jl schema reads.

An uncertainty is added alongside it whenever there is one to compute — that is, when the result
carries uncertainties at all (see [`calc_gradient_uncertainty`](@ref)). It costs about 2.4× the
gradient alone. When there is none the value is still added and `err` is left as it was, so
`val` carries one name that `err` does not; a table built from such results simply has no `σ_`
column for it.

```julia
result = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0)
with_gradient(result, :revelle, :fCO₂, :DIC; constant = :TA, relative = true).revelle
```
"""
function with_gradient(result::CarbonateResult, name::Symbol, numerator::Symbol,
                       denominator::Symbol; constant::Symbol, relative::Bool = false)
    haskey(result.val, name) && throw(ArgumentError(
        "this result already has a value called $name; choose another name"))

    value = relative ?
        calc_relative_gradient(result, numerator, denominator; constant = constant) :
        calc_gradient(result, numerator, denominator; constant = constant)

    values = merge(getfield(result, :val), NamedTuple{(name,)}((value,)))
    errors = _errors_with_gradient(result, name, numerator, denominator, constant, relative)

    return CarbonateResult(values, errors, result.input_keys, result.settings,
                           result.inputs, result.system, result.input_errors, result.Ks,
                           result.unit, result.source)
end

"`result.err` with the gradient's uncertainty added, where there is one to add."
function _errors_with_gradient(result::CarbonateResult, name::Symbol, numerator::Symbol,
                               denominator::Symbol, constant::Symbol, relative::Bool)
    errors = getfield(result, :err)
    isnothing(errors) && return errors

    uncertainty = calc_gradient_uncertainty(result, numerator, denominator;
                                            constant = constant, relative = relative)
    isnothing(uncertainty) && return errors

    return merge(errors, NamedTuple{(name,)}((uncertainty,)))
end
