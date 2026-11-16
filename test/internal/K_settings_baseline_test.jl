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
# and both must match the baseline. The one accepted difference: Millero 1979 no longer zeroes
# KB and the nutrient constants, which only ever mattered through a phosphate or silicate
# total that freshwater does not have (and where it did, gave PAlk = -PT). Those five
# constants are skipped for it; everything computed from them is still compared.

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

baseline_skips(name, field) =
    startswith(name, "constants_Millero 1979") && field in MILLERO1979_ZEROED

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
                for (how, call) in (("as written", kw), ("by hand", by_hand(kw)))
                    got = Dict(run_K_settings_case(kind, call))
                    wrong = [f for (f, v) in expected
                             if !baseline_skips(name, f) &&
                                !(haskey(got, f) && isapprox(got[f], v; rtol = 1e-12))]
                    isempty(wrong) || @info "$name ($how) differs from the baseline in" wrong
                    @test isempty(wrong)
                end
            end
        end
    end
end
