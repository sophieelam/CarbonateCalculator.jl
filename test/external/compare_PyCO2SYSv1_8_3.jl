# Agreement with PyCO2SYS v1.8.3, through `whole_system`. Baselines are made by
# generate_pyco2sys_baseline.py; the comparison itself is in agreement_common.jl.
using Test, CarbonateCalculator
isdefined(Main, :Agreement) || include(joinpath(@__DIR__, "agreement_common.jl"))
using .Agreement

@testset "PyCO2SYS v1.8.3: every option, S 33.1, 24 °C, 10 dbar (out: 12 °C, 1647 dbar)" begin
    run_agreement("PyCO2SYS v1.8.3 options",
                  joinpath(@__DIR__, "results", "compare_PyCO2SYS_v1_8_3.csv"),
                  CarbonateCalculator.whole_system)
end

@testset "PyCO2SYS v1.8.3: 0-35 °C, S 20-38, 0-11000 dbar (out: 2 °C, 4000 dbar)" begin
    run_agreement("PyCO2SYS v1.8.3 conditions",
                  joinpath(@__DIR__, "results", "compare_PyCO2SYS_v1_8_3_conditions.csv"),
                  CarbonateCalculator.whole_system)
end
