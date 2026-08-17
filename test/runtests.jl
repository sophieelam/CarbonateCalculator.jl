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

# The suite is in three blocks, each answering a different question and costing a different
# amount to run. `CC_TESTS` selects how far down the list to go:
#
#   quick     1 only        — does the package still agree with itself?
#   default   1 and 2       — and with other carbonate calculators? (the default; CI)
#   all       1, 2 and 3    — and with real-world measurements?
#
# Block 3 downloads published datasets, so it needs network access and a few hundred MB,
# which is why it is opt-in rather than skip-out.
#
#   CC_TESTS=quick julia --project=. test/runtests.jl
#   CC_TESTS=all CC_PLOTS=true julia --project=. test/runtests.jl
const TEST_LEVEL = get(ENV, "CC_TESTS", "default")
TEST_LEVEL in ("quick", "default", "all") ||
    error("CC_TESTS must be \"quick\", \"default\" or \"all\"; got \"$TEST_LEVEL\"")

const RUN_EXTERNAL_CALCULATORS = TEST_LEVEL in ("default", "all")
const RUN_FIELD_DATA = TEST_LEVEL == "all"


@testset "CarbonateCalculator.jl Full Test Suite" begin

    # --- 1. Internal consistency & regression -------------------------------------------
    # Everything the package can check about itself: unit behaviour of the speciation
    # functions, round trips, conservation across conditions, and uncertainty propagation.
    # No external reference data and no network, so this is the block to run on every edit.
    @testset "Internal consistency & regression" begin

        @testset "Carbon" begin
            include("internal/carbon_test.jl")
        end

        @testset "Boron" begin
            include("internal/boron_test.jl")
        end

        @testset "Boron Isotopes" begin
            include("internal/isotope_test.jl")
        end

        @testset "Input/Output Consistency" begin
            include("internal/input_output_test.jl")
        end

        @testset "Conditions" begin
            include("internal/conditions_test.jl")
        end

        @testset "pH Scales" begin
            include("internal/ph_scales_test.jl")
        end

        @testset "Paleo Proxies" begin
            include("internal/paleo_proxy_test.jl")
        end

        # Solve the same system from every valid pair of input parameters and require the
        # answers to agree — the package compared only against itself.
        @testset "Round-robin" begin
            include("internal/round_robin.jl")
            include("internal/cc_round_robin.jl")
        end

        @testset "Error Propagation" begin
            include("internal/errors_test.jl")
            include("internal/in_out_error_test.jl")
        end

    end

    # --- 2. Agreement with other carbonate calculators ----------------------------------
    # Saved PyCO2SYS and MATLAB CO2SYS outputs, compared parameter by parameter. These say
    # whether the package still implements the same chemistry as the reference tools, so
    # they are what matters when a numerical result changes. Roughly 55% of the runtime.
    if RUN_EXTERNAL_CALCULATORS
        @testset "External calculator agreement" begin
            include("external/compare_PyCO2SYSv1_8_3.jl")
            include("external/compare_MATLABv3_2_0.jl")
            include("external/compare_carbon_calculator.jl")
        end
    else
        @info "CC_TESTS=$TEST_LEVEL — skipping external calculator agreement."
    end

    # --- 3. Accuracy against real-world data --------------------------------------------
    # Published measurements where the carbonate system is over-determined, so a computed
    # parameter can be compared against the measured one. This is the only block that says
    # anything about accuracy rather than self-consistency.
    #
    # Set CC_PLOTS=true to also write comparison figures; off by default, since plotting
    # pulls in Plots and is a benchmarking aid rather than a check.
    if RUN_FIELD_DATA
        @testset "Real-world data accuracy" begin
            include("field_data/GLODAP_test.jl")
            include("field_data/CODAP_test.jl")
            include("field_data/SOCAT_test.jl")
            include("field_data/BGCArgo_test.jl")
        end
    else
        @info "CC_TESTS=$TEST_LEVEL — skipping real-world data accuracy " *
              "(set CC_TESTS=all; needs network access)."
    end

end