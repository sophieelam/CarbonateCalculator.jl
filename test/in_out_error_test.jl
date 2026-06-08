using Test
using ForwardDiff
using CarbonateCalculator

# --- SHARED TEST SETUP ---
# Define standard seawater inputs AND cold/deep output conditions
inputs = (
    TA = 2300.0,
    DIC = 2000.0,
    BT = 415.7,
    T_in = 25.0,
    P_in = 0.0,
    S_in = 35.0,
    T_out = 2.0,      # Deep water temperature
    P_out = 4000.0,   # High pressure at depth (dbar)
    S_out = 35.0,
    PT = 1.0,
    SiT = 15.0,
    δBT = 39.61,       # Adding explicit isotope for whole_system
    alphaB = 1.0272
)

# Define shared uncertainties
uncertainties = (
    TA = 2.0,
    DIC = 2.0,
    T_out = 0.5  
)

@testset "Carbonate Error Propagation Master Suite" begin

    @testset "Two-State: carbon_system" begin
        @info "Testing carbon_system (In/Out)..."
        result = carbon_system(; errors=uncertainties, K_method = "Lueker 2000", inputs...)
        
        @test hasproperty(result.err, :pHtot_in)
        @test hasproperty(result.err, :pHtot) 
        @test result.err.pHtot > 0.0
        
        # Error check: Out error should be larger due to T_out uncertainty
        @test result.err.pHtot > result.err.pHtot_in

        # Print the fancy table you like
        println("\n" * "═"^50)
        println("  TWO-STATE PROPAGATION RESULTS")
        println("  " * "─"^46)
        println("  Surface (In):")
        println("    pH Total:  $(round(result.val.pHtot_in, digits=4)) ± $(round(result.err.pHtot_in, digits=4))")
        println("\n  Deep Water (Out):")
        println("    pH Total:  $(round(result.val.pHtot, digits=4)) ± $(round(result.err.pHtot, digits=4))")
        println("═"^50 * "\n")
    end

    @testset "Isotopes: whole_system" begin
        @info "Testing whole_system isotopic errors..."
        # Add error for isotopic fractionation factor and delta
        iso_errs = merge(uncertainties, (δBT = 0.05, alphaB = 0.0001))
        
        result = whole_system(; errors=iso_errs, K_method="Lueker 2000", inputs...)
        
        @test hasproperty(result.err, :δBOH₄)
        @test result.err.δBOH₄ > 0.0
        @test result.err.pHtot > 0.0
    end

    @testset "Single-State: carbon_calculator" begin
        @info "Testing carbon_calculator..."
        result = carbon_calculator(; errors=(TA=2.0, DIC=2.0), inputs...)
        
        # Verify it's flat (No "_in" suffixes)
        @test hasproperty(result.err, :pHtot)
        @test !hasproperty(result.err, :pHtot_in)
        @test result.err.pHtot > 0.0
    end

    @testset "Boron Species: carbon_boron_calculator" begin
        @info "Testing carbon_boron_calculator..."
        result = carbon_boron_calculator(; errors=(BT=0.01,), inputs...)
        
        @test hasproperty(result.err, :BOH₃)
        @test result.err.BOH₃ > 0.0
    end

    @testset "Guardrails: MyAMI Intercepts" begin
        @info "Testing MyAMI Intercepts..."
        
        # Check carbon_system
        @test_logs (:warn, r"not supported") carbon_system(; 
            K_method="MyAMI", errors=(TA=2.0,), inputs...
        )

        # Check whole_system
        @test_logs (:warn, r"not supported") whole_system(; 
            K_method="MyAMI", errors=(TA=2.0,), inputs...
        )
    end
end