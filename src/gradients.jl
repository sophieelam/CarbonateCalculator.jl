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
    denominator === constant && throw(ArgumentError(
        "$denominator cannot be both varied and held constant"))

    scope = result.system
    base = _gradient_inputs(result, constant)
    start = getproperty(result.val, denominator)

    # The solver's own check, so an ill-posed pair is refused for the same reason and with the
    # same message it would get as input — including that pCO₂ and fCO₂ are one quantity rather
    # than two independent ones.
    _check_determinacy(merge(base, NamedTuple{(denominator,)}((start,))), scope;
                       require_two = scope === (:carbon,))

    solved_at(y) = getproperty(_solve_core(scope, merge(base, NamedTuple{(denominator,)}((y,)))),
                               numerator)

    return ForwardDiff.derivative(solved_at, start)
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

`err` is left untouched. A gradient is itself a derivative, so propagating input uncertainties
through it needs second-order AD; that is deliberately not attempted here.

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
    return CarbonateResult(values, result.err, result.input_keys, result.settings,
                           result.inputs, result.system, result.input_errors, result.Ks,
                           result.unit, result.is_collection_state)
end
