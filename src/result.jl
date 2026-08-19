# The data type returned by every calculation function.

"""
Keys carried in `CarbonateResult.settings`: everything that determines *how* the constants
were calculated, as opposed to the conditions or the chemistry.
"""
const SETTING_KEYS = (:K_method, :KSO4_method, :BT_method, :KF_method, :KNH3_method,
                      :Ca_method, :MyAMI_mode, :unit)

"""
Pick out the `SETTING_KEYS` from a call's inputs.

`NamedTuple{SETTING_KEYS}(inputs)` rather than filtering on `haskey`, so the key set is a
compile-time constant rather than something rebuilt on every call — worth several µs on a
path every calculation takes.

Raises on a missing key by design. Every caller builds its inputs from a signature that
declares all of these, so an absent one is a programming error rather than an input the user
left out, and should say so instead of quietly producing a smaller tuple.
"""
_settings(inputs::NamedTuple) = NamedTuple{SETTING_KEYS}(inputs)

"""
    CarbonateResult

The result of a carbonate system calculation.

Computed values are reached directly — `result.pHtot` — and forwarded to the underlying
`val` NamedTuple. `result.err` holds matching uncertainties when the calculation was run
with `errors=`, or `nothing`.

# Fields
- `val`: the computed state. **Numbers only** - every entry is a concentration, a condition
  or a derived quantity, so `result.pHtot` has a concrete type and arithmetic on it never
  meets a `nothing`. The two non-numeric pieces live in their own fields below.
- `Ks`: the equilibrium constants used.
- `unit`: the concentration unit the values are reported in.
- `err`: uncertainties matching `val`, or `nothing`.
- `input_keys`: which carbonate parameters the caller supplied, used by `show` to avoid
  reporting an input back as a result.
- `settings`: how the constants were calculated (see [`SETTING_KEYS`](@ref)).
- `inputs`: the original call, verbatim.
- `system`: which system was solved, `:carbon`, `:boron`, `:isotopes`, or `:whole`. Needed so a result can be
  re-solved without the caller having to say which function produced it.
- `input_errors`: the uncertainties this calculation was given, or `nothing`.

`settings`, `inputs`, `system` and `input_errors` exist so a result can be re-solved at other
conditions without the caller restating anything — see `at_collection_conditions`.
Keeping the *input* uncertainties, rather than only the propagated `err`, is what lets a
re-solve differentiate the whole chain from the original independent measurements. Propagating
`err` forward instead would double-count any input affecting both conditions, and would lose
the correlations between derived quantities.

`Ks` and `unit` are fields rather than entries in `val` so that `val` is numeric throughout,
which is what gives its fields concrete types without any parameterisation gymnastics.

The NamedTuples are held untyped rather than concretely parameterised so that
`ForwardDiff.Dual` values pass through unchanged.
"""
struct CarbonateResult{V<:NamedTuple, E<:Union{NamedTuple, Nothing}}
    val::V
    err::E
    input_keys::Vector{Symbol}
    settings::NamedTuple
    inputs::NamedTuple
    system::Tuple{Vararg{Symbol}}
    input_errors::Union{NamedTuple, Nothing}
    Ks::NamedTuple
    unit::String
end

"""
    CarbonateResult(raw, err, input_keys, settings, inputs, system, input_errors)

Build a result from a core's raw output, lifting `Ks` and `unit` out of the computed values.

The cores return the constants and the unit alongside the numbers, because both are needed to
finish the calculation; they become their own fields here so that `val` is numeric throughout.
Having this as a constructor rather than at each call site keeps the seven construction points
— four in the wrappers, three in `at_collection_conditions` — unaware of the split.
"""
function CarbonateResult(raw::NamedTuple, err, input_keys, settings, inputs, system,
                         input_errors)
    values, Ks, unit = _split_metadata(raw)
    return CarbonateResult(values, err, input_keys, settings, inputs, system, input_errors,
                           Ks, unit)
end

function Base.getproperty(res::CarbonateResult, s::Symbol)
    # 1. Standard fields
    if s in (:val, :err, :input_keys, :settings, :inputs, :system, :input_errors,
             :Ks, :unit)
        return getfield(res, s)
    end

    # Each result describes exactly one set of conditions, so a name means one thing.
    return getproperty(getfield(res, :val), s)
end

# Make tab-completion work for both the struct fields AND the math outputs
function Base.propertynames(res::CarbonateResult, private::Bool=false)
    return (fieldnames(CarbonateResult)..., keys(getfield(res, :val))...)
end

function Base.iterate(res::CarbonateResult, state...)
    return iterate(getfield(res, :val), state...)
end

function Base.keys(res::CarbonateResult)
    return keys(getfield(res, :val))
end

function Base.length(res::CarbonateResult)
    return length(getfield(res, :val))
end

# One result is one sample, so broadcast has to treat it as a scalar. Base's fallback is
# `broadcastable(x) = collect(x)`, which would reach the `iterate` above and splat a result
# into its ~38 values — making `f.(result, temperatures)` fail with a `DimensionMismatch`
# naming 38 axes instead of solving at each temperature.
Base.broadcastable(res::CarbonateResult) = Ref(res)