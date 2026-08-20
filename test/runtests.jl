using Test
using CarbonateCalculator

# Hand-checked reference values for the speciation functions, used by
# internal/{carbon,boron,isotope}_test.jl and nothing else. Included once here rather than
# per-file because it defines a module: three includes would evaluate `module CheckVals`
# three times and warn about replacing it.
include("internal/check_vals.jl")
using .CheckVals


# The suite is in three blocks, each answering a different question and costing a different
# amount to run. `CC_TESTS` selects how far down the list to go:
#
#   quick     1 only        — does the package still agree with itself?
#   default   1 and 2       — and with other carbonate calculators? (the default; CI)
#   all       1, 2 and 3    — and with real-world measurements?
#
# Block 3 downloads published datasets, so it needs network access and a few hundred MB,
# which is why it is opt-in rather than opt-out.
#
#   CC_TESTS=quick julia --project=. test/runtests.jl
#   CC_TESTS=all CC_PLOTS=true julia --project=. test/runtests.jl
#
# Block 1 is ~97% Julia compilation rather than computation: re-running a test file inside
# one session costs 1-4% of the first pass, and each distinct combination of supplied and
# omitted keywords is a fresh specialisation (~0.4s each). So `-O0` roughly halves the
# wall time and changes no result — worth it for the development loop:
#
#   CC_TESTS=quick julia -O0 --project=. test/runtests.jl     # 38s vs 72s
const TEST_LEVEL = get(ENV, "CC_TESTS", "default")
TEST_LEVEL in ("quick", "default", "all") ||
    error("CC_TESTS must be \"quick\", \"default\" or \"all\"; got \"$TEST_LEVEL\"")

const RUN_EXTERNAL_CALCULATORS = TEST_LEVEL in ("default", "all")
const RUN_FIELD_DATA = TEST_LEVEL == "all"

# Block 3's scaffolding is set up here rather than beside its testset below, because
# `using` and `const` are only legal at top level — inside the @testset they are a syntax
# error, and files brought in by `include` are evaluated in Main and cannot see a local.
if RUN_FIELD_DATA
    include("field_data/field_data_common.jl")
    using .FieldData

    # Plots is not a test dependency, so figures are skipped rather than fatal when
    # CC_PLOTS is set without it installed. The figure code lives in its own file because
    # it cannot even be *parsed* without Plots loaded — @layout is a macro, so an `if`
    # guard around it would not help.
    const FIELD_PLOTS = MAKE_PLOTS && optional_using(:Plots)
    FIELD_PLOTS && include("field_data/comparison_figure.jl")
end


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

        @testset "Input Validation" begin
            include("internal/validation_test.jl")
        end

        @testset "Method Selection" begin
            include("internal/method_selection_test.jl")
        end

        @testset "Core Consistency" begin
            include("internal/core_consistency_test.jl")
        end

        @testset "Solvers" begin
            include("internal/solver_test.jl")
        end

        @testset "Input/Output Consistency" begin
            include("internal/input_output_test.jl")
        end

        @testset "Tables.jl interop" begin
            include("internal/tables_test.jl")
        end

        @testset "Conditions" begin
            include("internal/conditions_test.jl")
        end

        @testset "Gradients" begin
            include("internal/gradient_test.jl")
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

        # Every `jldoctest` example in the docstrings, run and compared against the output it
        # claims. Last in the block because it is a documentation check rather than a
        # chemistry one.
        @testset "Doctests" begin
            include("internal/doctest_test.jl")
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
            include("field_data/BGCArgo_test.jl")
        end
    else
        @info "CC_TESTS=$TEST_LEVEL — skipping real-world data accuracy " *
              "(set CC_TESTS=all; needs network access)."
    end

end