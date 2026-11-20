# K_settings reproduces what the old flat keywords computed.
#
# test/refactor_baseline/K_settings_baseline.tsv holds every number the old code produced for
# the cases in K_settings_cases.jl, captured by capture_K_settings_baseline.jl from the source
# before K_settings existed. Each case is re-run here two ways:
#
#   as written    the old keywords, which the package now translates itself
#   by hand       the K_settings those keywords meant, spelled out below independently of the
#                 package's translation, so this checks it rather than repeating it
#
# and both must match the baseline, except where the package has since been corrected:
#
# 1. Millero 1979 no longer zeroes KB and the nutrient constants, which only ever mattered
#    through a phosphate or silicate total that freshwater does not have (and where it did,
#    gave PAlk = -PT). Those five constants are skipped for it; everything computed from them
#    is still compared.
# 2. Pressure corrections are now applied on the seawater scale (see `_to_total_mult` in
#    src/constants.jl), as in CO2SYS and PyCO2SYS. The old code corrected total-scale fits
#    (K1/K2 of Roy 1993, Lueker 2000, Sulpis 2020, SB 2020, Papadimitriou 2018; KB; KH2S; a
#    Clegg KNH3) directly on the total scale, which put pH about 1e-3 too high at 500 bar.
#    So at pressure, outside KGen (whose constants come from Kgen, unchanged), those five
#    constants and every solved quantity are not compared with the old numbers. The case is
#    still run both ways, and "as written" must still equal "by hand" to 1e-12: the
#    translation is what this test is for. Surface cases are unaffected (the correction is
#    exactly 1 at 0 bar) and are compared in full.
# 3. Four coefficient corrections found by comparison with CO2SYS v3 and PyCO2SYS (2026-11):
#    KH2S's salinity coefficient (0.0274, was 0.0247), WM13 KSO4 (a misplaced bracket),
#    Perez and Fraga KF (no mol/kg-H2O conversion) and aragonite's pressure ΔV (-45.96, was
#    -46.0). So KH2S is never compared; cases with H2ST, WM13 or Perez KF are compared only
#    for their conditions and totals; and KspA is not compared at pressure.
#
# Every exception is a skip, never a loosened tolerance, and "as written" = "by hand" is
# checked for every case regardless.

using Test, Logging, CarbonateCalculator

include(joinpath(@__DIR__, "..", "refactor_baseline", "K_settings_cases.jl"))

const K_BASELINE = joinpath(@__DIR__, "..", "refactor_baseline", "K_settings_baseline.tsv")

const OLD_METHOD_KEYWORDS = (:K_method, :KSO4_method, :BT_method, :KF_method, :KNH3_method,
                             :Ca_method, :MyAMI_mode)

"What each old K_method that switched companions meant, as K_settings."
const OLD_MEANING = Dict(
    "default" => (K1_K2_method = "KGen",),
    "Mehrbach 1973 A" => (K1_K2_method = "Mehrbach 1973", BT_method = "Culkin"),
    "Mehrbach 1973 B" => (K1_K2_method = "Mehrbach 1973", fH_method = "Peng 1987"),
    "Millero 1979" => (K1_K2_method = "Millero 1979", fH_method = "freshwater",
                       BT_method = "freshwater", pressure_method = "Millero 1983"),
)

"An old-style case rewritten with K_settings. The K_method meaning wins, as its override did."
function by_hand(kw::NamedTuple)
    flat = (; (k => v for (k, v) in pairs(kw) if k in OLD_METHOD_KEYWORDS && k !== :K_method)...)
    k = get(kw, :K_method, "default")
    settings = merge(flat, get(OLD_MEANING, k, (K1_K2_method = k,)))
    rest = (; (k => v for (k, v) in pairs(kw) if !(k in OLD_METHOD_KEYWORDS))...)
    return merge(rest, (K_settings = settings,))
end

const MILLERO1979_ZEROED = ("KB", "KP1", "KP2", "KP3", "KSi")

"Constants the corrections changed at pressure (differences 2 and 3 above)."
const PRESSURE_FIXED_CONSTANTS = ("K1", "K2", "KB", "KH2S", "KNH3", "KspA")
"Solved quantities that do not depend on the constants: conditions and totals."
const CONDITIONS_AND_TOTALS = ("temp_c", "sal", "pres_bar", "BT", "ST", "FT", "Ca", "Mg", "fH")

uses_kgen(kw) = get(kw, :K_method, "default") in ("default", "KGen")

"Was this case changed on purpose by the seawater-scale pressure correction (difference 2)?"
pressure_fixed(kw) = get(kw, :pres_bar, 0.0) > 0 && !uses_kgen(kw)

"Was every number in this case changed by a coefficient correction (difference 3)?"
corrected_throughout(kw) = get(kw, :H2ST, 0.0) > 0 ||
    (!uses_kgen(kw) && (get(kw, :KSO4_method, "") == "WM13" || get(kw, :KF_method, "") == "Perez"))

function baseline_skips(name, kind, kw, field)
    startswith(name, "constants_Millero 1979") && field in MILLERO1979_ZEROED && return true
    field == "KH2S" && return true                        # difference 3, at every condition
    corrected_throughout(kw) && return !(field in CONDITIONS_AND_TOTALS)
    pressure_fixed(kw) || return false       # KGen supplies its own KspA, unchanged
    return kind === :constants ? field in PRESSURE_FIXED_CONSTANTS :
                                 !(field in CONDITIONS_AND_TOTALS)
end

@testset "K_settings reproduces the old keywords" begin
    if !isfile(K_BASELINE)
        @warn "No K_settings baseline at $K_BASELINE; skipping. Capture it from the old " *
              "source with test/refactor_baseline/capture_K_settings_baseline.jl."
        @test_skip false
    else
        baseline = Dict{String, Dict{String, Float64}}()
        for line in Iterators.drop(eachline(K_BASELINE), 1)
            case, field, value = split(line, '\t')
            get!(baseline, case, Dict{String, Float64}())[field] = parse(Float64, value)
        end

        cases = K_settings_cases()
        @test Set(first.(cases)) == Set(keys(baseline))   # the case list has not drifted

        with_logger(NullLogger()) do
            for (name, (kind, kw)) in cases
                haskey(baseline, name) || continue
                expected = baseline[name]
                runs = Dict{String, Dict{String, Float64}}()
                for (how, call) in (("as written", kw), ("by hand", by_hand(kw)))
                    got = runs[how] = Dict(run_K_settings_case(kind, call))
                    wrong = [f for (f, v) in expected
                             if !baseline_skips(name, kind, kw, f) &&
                                !(haskey(got, f) && isapprox(got[f], v; rtol = 1e-12))]
                    isempty(wrong) || @info "$name ($how) differs from the baseline in" wrong
                    @test isempty(wrong)
                end
                # The translation itself, checked in full even where the baseline is not.
                a, b = runs["as written"], runs["by hand"]
                mismatch = [f for (f, v) in a
                            if !(haskey(b, f) && (isapprox(b[f], v; rtol = 1e-12) ||
                                                  (isnan(v) && isnan(b[f]))))]
                isempty(mismatch) || @info "$name: as written and by hand differ in" mismatch
                @test isempty(mismatch)
            end
        end
    end
end
