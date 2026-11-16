# Capture what the constants' settings produced *before* K_settings, so the change can be
# shown to reproduce them.
#
# Run this against the old source, then commit the TSV it writes:
#
#   git checkout main -- src                  # the code before K_settings
#   julia --project=. test/refactor_baseline/capture_K_settings_baseline.jl
#   git checkout HEAD -- src                  # back to the new code
#
# The internal test suite then re-runs every case on the new code (see
# test/internal/K_settings_baseline_test.jl). Uses the old API only, so it runs on both
# versions; on the new one it prints the deprecation warnings, which are silenced here.

using CarbonateCalculator, Logging

include(joinpath(@__DIR__, "K_settings_cases.jl"))

const OUT = joinpath(@__DIR__, isempty(ARGS) ? "K_settings_baseline.tsv" : ARGS[1])

open(OUT, "w") do io
    println(io, "case\tfield\tvalue")
    with_logger(NullLogger()) do
        for (name, (kind, kw)) in K_settings_cases()
            for (field, value) in run_K_settings_case(kind, kw)
                println(io, name, '\t', field, '\t', repr(value))
            end
        end
    end
end
println("wrote ", OUT)
