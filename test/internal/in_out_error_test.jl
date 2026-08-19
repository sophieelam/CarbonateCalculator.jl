using Test
using Logging
using ForwardDiff
using CarbonateCalculator

# --- SHARED TEST SETUP ---
# Conditions the sample was measured at.
inputs = (
    TA = 2300.0,
    DIC = 2000.0,
    BT = 415.7,
    temp_c = 25.0,
    pres_bar = 0.0,
    sal = 35.0,
    PT = 1.0,
    SiT = 15.0,
)

# The isotope arguments belong to whole_system only, so they are kept out of the shared tuple
# above. In it, carbon_system would absorb them into its kwargs sink and carry them into its
# output as fields that mean nothing for a carbon-only scope.
boron_inputs = merge(inputs, (
    δBT = 39.61,
    alphaB = 1.0272,
))

# Conditions the sample was collected at: cold and deep.
const TARGET_TEMP_C = 2.0
const TARGET_PRES_BAR = 400.0

# Uncertainties on the measurements themselves.
uncertainties = (
    TA = 2.0,
    DIC = 2.0,
)

# Propagation is reached by building a solver that names its uncertainties, then calling it
# positionally — the plain keyword API has no `errors` argument. These tests hold every
# parameter fixed and vary only the σ, so `varying` stays empty and the positional values
# are the uncertainties alone.
solve_with(σ::NamedTuple; scope=(:carbon,), settings...) =
    CarbonateSystem(scope...; varying_errors=keys(σ), settings...)(values(σ)...)

@testset "Carbonate Error Propagation Master Suite" begin

    @testset "Two-stage: carbon_system then recalculate" begin
        @info "Testing carbon_system across measured and collection conditions..."
        measured = solve_with(uncertainties; K_method="Lueker 2000", inputs...)
        @test measured.err.pHtot > 0.0

        # Measurement uncertainties are carried on the result, so they are not restated -
        # the re-solve differentiates the whole chain from the original measurements.
        at_target = recalculate_at_target_conditions(measured;
            temp_c=TARGET_TEMP_C, pres_bar=TARGET_PRES_BAR)
        @test at_target.err.pHtot > 0.0
        @test isfinite(at_target.err.pHtot)

        # Adding uncertainty in the collection temperature can only widen the result.
        with_target_err = recalculate_at_target_conditions(measured;
            temp_c=TARGET_TEMP_C, pres_bar=TARGET_PRES_BAR, errors=(temp_c=0.5,))
        @test with_target_err.err.pHtot > at_target.err.pHtot

        # A result with no uncertainties attached stays that way.
        no_errors = carbon_system(; K_method="Lueker 2000", inputs...)
        @test isnothing(recalculate_at_target_conditions(no_errors; temp_c=TARGET_TEMP_C).err)

        println("\n" * "═"^50)
        println("  TWO-STAGE PROPAGATION RESULTS")
        println("  " * "─"^46)
        println("  Measured (25 °C, 0 bar):")
        println("    pH Total:  $(round(measured.val.pHtot, digits=4)) ± $(round(measured.err.pHtot, digits=4))")
        println("\n  Collected ($(TARGET_TEMP_C) °C, $(TARGET_PRES_BAR) bar):")
        println("    pH Total:  $(round(at_target.val.pHtot, digits=4)) ± $(round(at_target.err.pHtot, digits=4))")
        println("    ...with ±0.5 °C on collection temperature: ± $(round(with_target_err.err.pHtot, digits=4))")
        println("═"^50 * "\n")
    end

    @testset "Isotopes: whole_system" begin
        @info "Testing whole_system isotopic errors..."
        # Add error for isotopic fractionation factor and delta
        iso_errs = merge(uncertainties, (δBT = 0.05, alphaB = 0.0001))

        result = solve_with(iso_errs; scope=(:carbon, :boron, :isotopes),
                            K_method="Lueker 2000", boron_inputs...)

        @test hasproperty(result.err, :δBOH₄)
        @test result.err.δBOH₄ > 0.0
        @test result.err.pHtot > 0.0
    end

    @testset "Single-State: carbon_calculator" begin
        @info "Testing carbon_calculator..."
        result = solve_with((TA=2.0, DIC=2.0); inputs...)

        # One result describes one set of conditions, so there are no _in/_out names.
        @test hasproperty(result.err, :pHtot)
        @test !hasproperty(result.err, :pHtot_in)
        @test result.err.pHtot > 0.0
    end

    @testset "Boron Species: carbon_boron_calculator" begin
        @info "Testing carbon_boron_calculator..."
        result = solve_with((BT=0.01,); scope=(:carbon, :boron, :isotopes), boron_inputs...)
        
        @test hasproperty(result.err, :BOH₃)
        @test result.err.BOH₃ > 0.0
    end

    @testset "MyAMI supports error propagation" begin
        @info "Testing error propagation through the MyAMI path..."
        # Kgen.jl is pure Julia and Real-typed, so ForwardDiff differentiates straight through
        # it. Anything on this path that is not AD-transparent discards the requested errors,
        # which shows up as uncertainties of exactly zero rather than as a failure.

        for (scope, args) in (((:carbon,), inputs),
                              ((:carbon, :boron, :isotopes), boron_inputs))
            result = solve_with((TA=2.0,); scope=scope, K_method="MyAMI", args...)
            @test hasproperty(result, :err)
            @test result.err.pHtot > 0.0
            @test isfinite(result.err.pHtot)
        end

        # No warning should be emitted any more.
        @test_logs min_level=Logging.Warn solve_with((TA=2.0,);
            K_method="MyAMI", inputs...
        )
    end

    @testset "Guardrails: MyAMI_mode" begin
        @info "Testing MyAMI_mode validation..."
        # Only the polynomial approximation is implemented in Kgen.jl; the full MyAMI
        # model is Python-only. Asking for it must fail loudly rather than silently
        # returning approximate values.
        @test_throws ArgumentError carbon_system(;
            K_method="MyAMI", MyAMI_mode="calculate", inputs...
        )
        @test carbon_system(; K_method="MyAMI", MyAMI_mode="approximate", inputs...).pHtot ≈
              carbon_system(; K_method="MyAMI", inputs...).pHtot
    end
end