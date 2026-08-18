# The data type returned by every calculation function.

"""
Keys carried in `CarbonateResult.settings`: everything that determines *how* the constants
were calculated, as opposed to the conditions or the chemistry.
"""
const SETTING_KEYS = (:K_method, :KSO4_method, :BT_method, :KF_method, :KNH3_method,
                      :Ca_method, :MyAMI_mode, :unit, :scale)

"Pick out the `SETTING_KEYS` that are actually present in a call's inputs."
function _settings(inputs::NamedTuple)
    present = Tuple(k for k in SETTING_KEYS if haskey(inputs, k))
    return NamedTuple{present}(Tuple(inputs[k] for k in present))
end

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
- `val`: the computed state.
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

The NamedTuples are held untyped rather than concretely parameterised so that
`ForwardDiff.Dual` values pass through unchanged.
"""
struct CarbonateResult
    val::NamedTuple
    err::Union{NamedTuple, Nothing}
    input_keys::Vector{Symbol}
    settings::NamedTuple
    inputs::NamedTuple
    system::Symbol
    input_errors::Union{NamedTuple, Nothing}
end

function Base.getproperty(res::CarbonateResult, s::Symbol)
    # 1. Standard fields
    if s in (:val, :err, :input_keys, :settings, :inputs, :system, :input_errors)
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