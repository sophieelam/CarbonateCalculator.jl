# Agreement with CO2SYS-MATLAB v3.2.0, through `whole_system`. The baseline is made by
# compare_MATLABv3_2_0.m; the comparison itself is in agreement_common.jl.
#
# PyCO2SYS v1.8.3 reproduces every row of this baseline to 1.4e-7 in pH and 4e-7 relative
# in DIC, CO3 and fCO2, on all four pH scales and every option, so the two references
# confirm each other: a disagreement here is a difference in this package.
using Test, CarbonateCalculator
isdefined(Main, :Agreement) || include(joinpath(@__DIR__, "agreement_common.jl"))
using .Agreement

@testset "MATLAB CO2SYS v3.2.0: every option, S 33.1, 24 °C, 1 dbar (out: 12 °C, 1647 dbar)" begin
    run_agreement("MATLAB CO2SYS v3.2.0",
                  joinpath(@__DIR__, "results", "compare_MATLABv3_2_0.csv"),
                  CarbonateCalculator.whole_system)
end
