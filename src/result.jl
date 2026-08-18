# The data type returned by every calculation function.

"""
Keys carried in `CarbonateResult.settings`: everything that determines *how* the constants
were calculated, as opposed to the conditions or the chemistry.
"""
const SETTING_KEYS = (:K_method, :KSO4_method, :BT_method, :KF_method, :KNH3_method,
                      :Ca_method, :MyAMI_mode, :unit)

"""
Pick out the `SETTING_KEYS` from a call's inputs.

`NamedTuple{SETTING_KEYS}(inputs)` rather than filtering on `haskey`: the key set is then a
compile-time constant instead of a value computed per call, which took this from 4.1 µs to
negligible. Every caller — the four wrappers and `recalculate_at_target_conditions` — builds
its inputs from a signature that declares all of these, so a missing key means a genuine
programming error and should say so rather than silently produce a smaller tuple.
"""
_settings(inputs::NamedTuple) = NamedTuple{SETTING_KEYS}(inputs)

"Arguments that were removed, and what to use instead."
const RETIRED_ARGUMENTS = (
    T_in  = "temp_c",
    S_in  = "sal",
    P_in  = "pres_bar",
    T_out = "recalculate_at_target_conditions(result; temp_c=...)",
    S_out = "nothing - a sample's salinity does not change between collection and measurement",
    P_out = "recalculate_at_target_conditions(result; pres_bar=...)",
    # Accepted and documented for years, threaded through every signature, and never once
    # applied - `carbon_system(..., pdict=Dict(:temp_c=>2.0))` returned the same answer as
    # omitting it. Julia splatting does the job natively and composes better, so the
    # argument is gone rather than implemented.
    pdict = "splatting: carbon_system(; TA=2300.0, DIC=2000.0, your_parameters...)",
    # Chose between one K_method for a whole array and one per sample - but only inside
    # K_calculator's array branch, which the public entry points could never reach because
    # they reject arrays first. Both are gone; process many samples with a solver instead.
    K_mode = "carbon_solver, which builds a scalar solver you can broadcast",
    # A generic pH plus `scale` said exactly what the scale-specific arguments already say,
    # at the cost of a permanently-`nothing` `pH` output field and a setting that existed
    # only to disambiguate it. Nothing used it.
    pH = "pHtot, pHsws, pHfree or pHNBS - name the scale directly",
    scale = "nothing - name the scale in the argument, e.g. pHsws=8.0",
)

"""
Reject arguments that no longer exist.

Every calculation function ends in `kwargs...`, so an unrecognised keyword is silently
absorbed rather than raising. Without this check, a call still written against the old API
would quietly compute at the default 25 °C and return a plausible wrong number.
"""
function _reject_retired_arguments(kwargs)
    found = [k for k in keys(kwargs) if haskey(RETIRED_ARGUMENTS, k)]
    isempty(found) && return nothing

    lines = ["$k is no longer accepted; use $(RETIRED_ARGUMENTS[k])" for k in found]

    # The condition rename needs a word of explanation that the others do not.
    conditions = [k for k in found if k in (:T_in, :S_in, :P_in, :T_out, :S_out, :P_out)]
    note = isempty(conditions) ? "" :
        "\nConditions are now named temp_c, sal and pres_bar, and a single call describes " *
        "one set of conditions."

    throw(ArgumentError("retired argument(s):\n  " * join(lines, "\n  ") * note))
end

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
- `system`: which system was solved, `:carbon` or `:whole`. Needed so a result can be
  re-solved without the caller having to say which function produced it.
- `input_errors`: the uncertainties this calculation was given, or `nothing`.

`settings`, `inputs`, `system` and `input_errors` exist so a result can be re-solved at
other conditions without the caller restating anything — see
`recalculate_at_target_conditions`. Keeping the *input* uncertainties, rather than only the
propagated `err`, is what lets a re-solve differentiate the whole chain from the original
independent measurements; propagating `err` forward instead would double-count any input
that affects both conditions, and would lose the correlations between derived quantities.

`val` used to carry `Ks`, `unit`, and a `pH` field that was `nothing` on every call - the
generic pH slot the core cleared after mapping it onto `pHtot`. Splitting the two metadata
entries out and retiring the generic `pH` leaves `val` numeric throughout, which is what makes
its field types concrete without any parameterisation gymnastics.

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
— four in the wrappers, three in `recalculate_at_target_conditions` — unaware of the split.
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

# Allow iteration so tools like ForwardDiff can treat the result like a tuple
function Base.iterate(res::CarbonateResult, state...)
    return iterate(getfield(res, :val), state...)
end

function Base.keys(res::CarbonateResult)
    return keys(getfield(res, :val))
end

function Base.length(res::CarbonateResult)
    return length(getfield(res, :val))
end