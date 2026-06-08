using Test
# If running this outside your main package environment, uncomment the line below:
using CarbonateCalculator 

@testset "Automatic Differentiation Error Propagation" begin
    @info "Testing Error Propagation via ForwardDiff..."

    # 1. Define standard seawater inputs
    base_inputs = (
        TA = 2300.0,
        DIC = 2000.0,
        T_in = 25.0,
        S_in = 35.0,
        PT = 1.0,
        SiT = 15.0,
        unit = "umol",
        scale = "total"
    )

    # 2. Define realistic measurement uncertainties
    # e.g., ±2.0 μmol/kg for TA and DIC, ±0.1 for Silicate
    uncertainties = (
        TA = 2.0,
        DIC = 2.0,
        SiT = 0.1
    )

    # 3. Call the error propagator
    # Note: We pass the `whole_system` function itself as the first argument
    result = propagate_errors(whole_system, inputs=base_inputs, errors=uncertainties)

    vals = result.val
    errs = result.err

    # --- SANITY CHECKS ---
    
    # Ensure the returned structure contains the expected keys
    @test hasproperty(errs, :pHtot)
    @test hasproperty(errs, :pCO₂)
    @test hasproperty(errs, :ΩC)
    @test hasproperty(errs, :ΩA)

    # Errors must be strictly positive real numbers (no NaNs, no negatives)
    @test errs.pHtot > 0.0
    @test errs.pCO₂ > 0.0
    @test errs.ΩC > 0.0
    @test !isnan(errs.pHtot)

    # --- MATH CHECKS ---
    
    # The propagated error of the inputs themselves should exactly match 
    # the input uncertainties (since ∂TA/∂TA = 1)
    @test isapprox(errs.TA, 2.0, atol=1e-8)
    @test isapprox(errs.DIC, 2.0, atol=1e-8)
    @test isapprox(errs.SiT, 0.1, atol=1e-8)

    # Variables independent of the inputs with errors should have 0 uncertainty
    # (Since temperature has no input uncertainty defined here, it shouldn't have an error)
    if hasproperty(errs, :T_in)
        @test isapprox(errs.T_in, 0.0, atol=1e-12)
    end

    # Print a nice summary for the terminal
    println("\n" * "═"^40)
    println("  PROPAGATION RESULTS")
    println("  " * "─"^36)
    println("  pH Total:  $(round(vals.pHtot, digits=4)) ± $(round(errs.pHtot, digits=4))")
    println("  pCO₂:      $(round(vals.pCO₂, digits=1)) ± $(round(errs.pCO₂, digits=1)) μatm")
    println("  Ω Calcite: $(round(vals.ΩC, digits=2)) ± $(round(errs.ΩC, digits=2))")
    println("═"^40 * "\n")
end