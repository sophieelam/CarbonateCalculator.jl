using Test
using CarbonateCalculator
include("check_vals.jl")
using .CheckVals

import .CarbonateCalculator.Carbon: CO₂_from_pH_DIC, H_from_HCO₃_CO₃, H_from_HCO₃_TA, 
                           pH_from_HCO₃_DIC, H_from_CO₃_TA, H_from_CO₃_DIC, 
                           pH_from_TA_DIC, calc_CO₂, calc_CO₃, calc_HCO₃, calc_TA, 
                           fCO₂_to_CO₂, CO₂_to_fCO₂, fCO₂_to_pCO₂, pCO₂_to_fCO₂,
                           DIC_from_CO₂_pH, H_from_CO₂_HCO₃, H_from_CO₂_CO₃,
                           pH_from_CO₂_TA, H_from_CO₂_DIC, DIC_from_pH_HCO₃,
                           DIC_from_pH_CO₃, DIC_from_pH_TA
import .CarbonateCalculator.Boron: H_from_BT_BOH3, H_from_BT_BOH4, BT_from_pH_BOH3, 
                          BT_from_pH_BOH4, calc_BOH3, calc_BOH4, calc_chiB,
                          B_calculator
import .CarbonateCalculator.Isotopes: get_alphaB, calc_ABT, ABOH3_from_H_ABT, ABOH4_from_H_ABT,
                            A11_to_δ11, δ11_to_A11


@testset "CarbonateCalculator.jl Full Test Suite" begin

    @testset "Carbon" begin
        include("carbon_test.jl")
    end

    @testset "Boron" begin
        include("boron_test.jl")
    end

    @testset "Boron Isotopes" begin
        include("isotope_test.jl")
    end

    @testset "Input/Output Consistency" begin
        include("input_output_test.jl")
    end

    @testset "Paleo Proxies" begin
        include("paleo_proxy_test.jl")
    end

    @testset "CO2SYS Replication & Consistency" begin
        include("CO2SYS_tests/round_robin.jl")
        include("CO2SYS_tests/cc_round_robin.jl")
        include("CO2SYS_tests/compare_PyCO2SYSv1_8_3.jl")
        include("CO2SYS_tests/compare_MATLABv3_2_0.jl")
        include("CO2SYS_tests/compare_carbon_calculator.jl")
    end

        @testset "Error Propagation" begin
        include("errors_test.jl")
        include("in_out_error_test.jl")
    end

end