# Verification step 2: the two-stage API must reproduce the pre-refactor `_out` values,
# except where the old recursion is known to have dropped an input.
#
#   julia --project=. test/refactor_baseline/check_against_baseline.jl
#
# Reads `baseline_two_condition.tsv`, re-runs each case as stage 1 + stage 2, and compares.
# Exits non-zero if anything differs that is not an accounted-for bug fix.

using CarbonateCalculator
using Printf

const BASELINE = joinpath(@__DIR__, "baseline_two_condition.tsv")
const TOL = 1e-9
const δBT_SW = 39.61

# The baseline was captured before Helpers.calc_ST/calc_BT were deleted in favour of the
# Constants versions, so it embeds the old total sulfate (which divided by 90.062 rather
# than the 96.062 molar mass of SO₄²⁻, 6.7% high) and the old total boron (0.000416 rather
# than Uppstrom's 0.0004157). It cannot be regenerated - the two-condition code is gone.
#
# Pinning those totals here keeps the comparison about what it is meant to test: whether the
# two-stage API reproduces the old *conditions* behaviour. Without them every case differs
# by the totals fix and a genuine regression would be lost in the noise.
# Only the MyAMI branch used the Helpers formulas; the analytic K_methods already used the
# Constants ones, so their baseline values need no pinning. That split is itself the bug.
const OLD_ST = 0.14 * 35.0 / 1.80655 / 90.062
const OLD_FT = 6.7e-5 * 35.0 / 1.80655 / 18.9984
const OLD_BT_UMOL = (0.000416 / 35.0) * 35.0 * 1e6

uses_myami(kw) = !haskey(kw, :K_method) || kw.K_method == "MyAMI"

"The totals the old code would have used for this case, BT in the case's own unit."
function old_totals(kw)
    uses_myami(kw) || return NamedTuple()
    return (ST = OLD_ST,
            FT = OLD_FT,
            BT = OLD_BT_UMOL / (get(kw, :unit, "umol") == "mmol" ? 1e3 : 1.0))
end

# Cases whose target-condition state the old code got wrong. Each entry says why.
const KNOWN_BUGS = Dict(
    "csys_H2S_NH4"       => "H2ST/NH4T were never forwarded, so target-condition alkalinity omitted them",
    "csys_lueker_khoo"   => "KSO4_method was hardcoded to \"default\" in the recursion",
    "wsys_dBT_nondefault" => "δBT was never forwarded, so the target state re-derived the 39.61 default",
)

# Cases that no longer exist: changing salinity between collection and measurement is not
# physical, and is now rejected rather than silently rescaling BT/ST/FT.
const RETIRED_CASES = ("csys_S_28to38", "wsys_S_28to38")

# case => (stage-1 kwargs, target temp_c, target pres_bar)
const CASES = Dict(
    "csys_T_20to2"        => ((TA=2300.0, DIC=2000.0, temp_c=20.0), 2.0, nothing),
    "csys_P_0to400"       => ((TA=2300.0, DIC=2000.0, temp_c=20.0, pres_bar=0.0), nothing, 400.0),
    "csys_TP_both"        => ((TA=2300.0, DIC=2000.0, temp_c=20.0, pres_bar=0.0), 2.0, 400.0),
    "csys_pHtot_TA"       => ((pHtot=8.1, TA=2300.0, temp_c=25.0), 5.0, nothing),
    "csys_pCO2_TA"        => ((pCO₂=400.0, TA=2300.0, temp_c=25.0), 5.0, nothing),
    "csys_nutrients"      => ((TA=2300.0, DIC=2000.0, PT=1.0, SiT=15.0, temp_c=25.0), 5.0, nothing),
    "csys_H2S_NH4"        => ((TA=2300.0, DIC=2000.0, H2ST=5.0, NH4T=3.0, temp_c=25.0), 5.0, nothing),
    "csys_lueker_khoo"    => ((TA=2300.0, DIC=2000.0, temp_c=20.0, K_method="Lueker 2000", KSO4_method="Khoo"), 2.0, nothing),
    "csys_lueker_default" => ((TA=2300.0, DIC=2000.0, temp_c=20.0, K_method="Lueker 2000"), 2.0, nothing),
    "csys_unit_mmol"      => ((TA=2.3, DIC=2.0, unit="mmol", temp_c=20.0), 2.0, nothing),
    "csys_paleo_CaMg"     => ((TA=2300.0, DIC=2000.0, Ca=0.02, Mg=0.03, temp_c=20.0), 2.0, nothing),
    "wsys_dBT_T_20to2"    => ((TA=2300.0, DIC=2000.0, δBT=δBT_SW, temp_c=20.0), 2.0, nothing),
    "wsys_pHtot_dBT"      => ((pHtot=8.1, TA=2300.0, δBT=δBT_SW, temp_c=20.0), 2.0, nothing),
    "wsys_TP_both"        => ((TA=2300.0, DIC=2000.0, δBT=δBT_SW, temp_c=20.0, pres_bar=0.0), 2.0, 400.0),
    "wsys_alphaB"         => ((TA=2300.0, DIC=2000.0, δBT=δBT_SW, alphaB=1.0272, temp_c=20.0), 2.0, nothing),
    "wsys_dBT_nondefault" => ((TA=2300.0, DIC=2000.0, δBT=38.0, temp_c=20.0), 2.0, nothing),
)

is_whole_system(case) = startswith(case, "wsys")

"Read the baseline into case => Dict(key => value)."
function read_baseline()
    out = Dict{String,Dict{String,Float64}}()
    for (i, line) in enumerate(eachline(BASELINE))
        i == 1 && continue
        case, key, value = split(line, '\t')
        get!(out, case, Dict{String,Float64}())[key] = parse(Float64, value)
    end
    return out
end

baseline = read_baseline()
unexplained = 0
bugfixes = 0
skipped = 0

for case in sort(collect(keys(baseline)))
    if case in RETIRED_CASES
        @printf("%-22s skipped (salinity change no longer supported)\n", case)
        global skipped += 1
        continue
    end
    haskey(CASES, case) || continue

    kw, tc, pb = CASES[case]
    f = is_whole_system(case) ? whole_system : carbon_system
    pinned = merge(kw, old_totals(kw))
    new = at_collection_conditions(f(; pinned...); temp_c=tc, pres_bar=pb).val

    worst_key, worst_rel = "", 0.0
    for (key, old) in baseline[case]
        endswith(key, "_out") || continue
        base = Symbol(chopsuffix(key, "_out"))
        haskey(new, base) || continue
        n = new[base]
        (n isa Number && isfinite(n)) || continue
        rel = abs(old - n) / max(abs(old), 1e-30)
        rel > worst_rel && ((worst_key, worst_rel) = (key, rel))
    end

    if worst_rel <= TOL
        @printf("%-22s match (worst rel=%.1e)\n", case, worst_rel)
    elseif haskey(KNOWN_BUGS, case)
        @printf("%-22s DIFFERS as expected — %s\n", case, KNOWN_BUGS[case])
        @printf("%22s   worst %s rel=%.3e\n", "", worst_key, worst_rel)
        global bugfixes += 1
    else
        @printf("%-22s DIFFERS UNEXPECTEDLY — worst %s rel=%.3e\n", case, worst_key, worst_rel)
        global unexplained += 1
    end
end

println()
println("matched            : ", length(baseline) - bugfixes - unexplained - skipped)
println("expected bug fixes : ", bugfixes)
println("skipped (retired)  : ", skipped)
println("unexplained        : ", unexplained)
exit(unexplained == 0 ? 0 : 1)
