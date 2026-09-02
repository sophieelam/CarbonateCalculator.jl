# src/errors.jl
using ForwardDiff, LinearAlgebra

"""
Substitute the perturbed values back into the inputs, keeping the result concretely typed.

`Val` carries the key set as a type parameter, so the NamedTuple is built from a compile-time
key list and a fixed-length tuple. `Tuple(values)` would not do: the length of a `Vector` is not
known to the compiler, so the result would be a `Vararg` tuple and the types would be lost
again, entering the solver as `Any` on every Dual pass.
"""
_perturbed(inputs::NamedTuple, ::Val{names}, values) where {names} =
    merge(inputs, NamedTuple{names}(ntuple(i -> values[i], Val(length(names)))))

"""
The numeric entries of a solved state, as the vector ForwardDiff differentiates.

`@generated` so the numeric fields are chosen from the NamedTuple's *type*: a solved state also
carries the equilibrium constants as a nested NamedTuple and the unit as a `String`, and
filtering those out with a runtime `v isa Number` would build a `Vector{Any}` on every
differentiated pass.
"""
@generated function _numeric_values(state::NamedTuple{names, types}) where {names, types}
    numeric = [:(state.$name) for (name, T) in zip(names, types.parameters) if T <: Number]
    return :([$(numeric...)])
end

"Names matching [`_numeric_values`](@ref), so the Jacobian's rows can be labelled."
@generated function _numeric_names(::NamedTuple{names, types}) where {names, types}
    return :($(Tuple(name for (name, T) in zip(names, types.parameters) if T <: Number)))
end

"""
Drop the derivative parts of a solved state, recovering the values the Jacobian pass computed.

A Dual carries its primal through arithmetic untouched, so these are the same numbers a separate
`Float64` evaluation would produce, and taking them from the Jacobian pass saves one solve per
row.
"""
_strip_partials(x::ForwardDiff.Dual) = ForwardDiff.value(x)
_strip_partials(x::NamedTuple) = map(_strip_partials, x)   # recurses into the constants
_strip_partials(x) = x

"""
The reported values, taken from the differentiated pass where that is possible.

`propagate_errors` accepts any function of the inputs. The solvers in this package return a raw
solved state, which is a NamedTuple all the way down and so can have its derivative parts
stripped. A public entry point returns a `CarbonateResult`, and rebuilding one of those field by
field — `input_keys`, `inputs`, `settings` and all — is more surface area than the solve it
would save, so that case is simply re-solved.
"""
_recover_values(state::NamedTuple, target_func, inputs) = _strip_partials(state)
_recover_values(::Any, target_func, inputs) = target_func(; inputs...)

"""
    propagate_errors(target_func; inputs::NamedTuple, errors::NamedTuple)

Propagate uncertainties through a carbonate chemistry function by first-order Taylor expansion,
with the derivatives from automatic differentiation.

`target_func` is called with `inputs` as keyword arguments; `errors` gives a σ for whichever of
them are uncertain. Returns `(val, err)` — the solved state, and a matching σ for every numeric
quantity in it. Because the derivatives are taken rather than written down, every derived
quantity gets an uncertainty without anyone deriving an error formula for it.

!!! warning
    Assumes the named inputs are uncorrelated. See below.

# Extended help

## The model, and its one assumption

σ_out = √( Σᵢ (∂out/∂inᵢ · σᵢ)² ), which **assumes the named inputs are uncorrelated**. There is
no covariance term, so uncertainties that share a common cause — the same instrument, a shared
calibration offset — are not represented, and would need a full covariance treatment.

That assumption is about the *named inputs*, not about the routes through the calculation. An
input reaching the answer several ways is handled exactly, because it is differentiated once:
salinity enters every equilibrium constant, and TA and DIC enter both conditions of a two-stage
calculation, with the correlations that implies preserved by construction. This is why
[`at_collection_conditions`](@ref) differentiates the whole chain from the original measurements
rather than propagating an intermediate σ forward — the latter would double-count anything
affecting both conditions.

First-order also means the result is a local linearisation. For a σ large enough that the
system is meaningfully non-linear across it, it is an approximation.
"""
function propagate_errors(target_func; inputs::NamedTuple, errors::NamedTuple)
    error_keys = keys(errors)

    # Must stay above the `Float64[...]` conversion below. Underneath it this could never
    # fire, because the conversion meets the `nothing` first and raises
    # `MethodError: Cannot convert an object of type Nothing to an object of type Float64`,
    # which names neither the parameter nor the reason.
    for key in error_keys
        isnothing(getproperty(inputs, key)) && throw(ArgumentError(
            "an uncertainty was given for $key, but no value for $key was supplied.\n" *
            "Uncertainties propagate from the inputs you measured, so $key needs a value " *
            "before it can have an error."))
    end

    # Ensure we are working with Floats for the AD process
    base_values = Float64[getproperty(inputs, k) for k in error_keys]
    uncertainties = Float64[getproperty(errors, k) for k in error_keys]

    # The last state the Jacobian pass solved, kept so its values can be reused rather than
    # recomputed. Every pass produces the same primal, so which one it ends up holding does not
    # matter — including under ForwardDiff's chunked mode, where there is more than one.
    solved = Ref{Any}(nothing)

    error_names = Val(error_keys)
    function wrapped_math(x_vec)
        state = target_func(; _perturbed(inputs, error_names, x_vec)...)
        solved[] = state
        return _numeric_values(state)
    end

    # The matrix of partial derivatives, one row per output, one column per uncertain input.
    jac = ForwardDiff.jacobian(wrapped_math, base_values)

    # σ_out = sqrt( Σ (∂out/∂in * σ_in)^2 ), summed across the input contributions.
    out_variances = (jac .^ 2) * (uncertainties .^ 2)
    out_errors_vec = sqrt.(sum(out_variances, dims=2))

    baseline_state = _recover_values(solved[], target_func, inputs)
    err_tuple = NamedTuple{_numeric_names(baseline_state)}(vec(out_errors_vec))

    return (val = baseline_state, err = err_tuple)
end
