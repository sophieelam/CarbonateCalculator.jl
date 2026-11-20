# The same comparisons through `carbon_system`, the carbonate-only preset, which must give
# the same answers as `whole_system` for the carbonate system. See agreement_common.jl.
using Test, CarbonateCalculator
isdefined(Main, :Agreement) || include(joinpath(@__DIR__, "agreement_common.jl"))
using .Agreement

@testset "carbon_system" begin
    @testset "against MATLAB CO2SYS v3.2.0" begin
        run_agreement("carbon_system vs MATLAB CO2SYS v3.2.0",
                      joinpath(@__DIR__, "results", "compare_MATLABv3_2_0.csv"),
                      CarbonateCalculator.carbon_system)
    end
    @testset "against PyCO2SYS v1.8.3, all conditions" begin
        run_agreement("carbon_system vs PyCO2SYS v1.8.3 conditions",
                      joinpath(@__DIR__, "results", "compare_PyCO2SYS_v1_8_3_conditions.csv"),
                      CarbonateCalculator.carbon_system)
    end
end
