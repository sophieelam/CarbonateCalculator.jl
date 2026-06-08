using Test

# 1. Define the baseline environmental conditions 
ENV_KWARGS = (
    S_in = 33.0, 
    T_in = 22.0, 
    P_in = 1234.0, 
    Si_in = 10.0,
    PO4_in = 1.0,
    NH3_in = 2.0,
    H2S_in = 3.0,
    K_method = "Lueker 2000"
)

@testset "Round-Robin Internal Consistency" begin

    # 2. Solve the base system using TA and DIC
    base_TA = 2300.0
    base_DIC = 2100.0
    
    # Calculate the baseline truth
    baseline = whole_system(TA=base_TA, DIC=base_DIC; ENV_KWARGS...)

    # 3. Map out the parameter names EXACTLY as Calculator.jl returns them
    calculated_vars = Dict(
        :TA => baseline.TA,
        :DIC => baseline.DIC,
        :pHtot => baseline.pHtot,   # Changed from :pH
        :pCO₂ => baseline.pCO₂,     # Note the subscript '2' - check Calculator.jl exports!
        :CO₃ => baseline.CO₃,
        :HCO₃ => baseline.HCO₃,
        :CO₂ => baseline.CO₂
    )

    all_pars = collect(keys(calculated_vars))

    # 4. Identify Invalid Pairs
    # Added pCO2 & CO2 because they are chemically redundant
    invalid_pairs = Set([
        Set([:pCO₂, :CO₂]), 
        Set([:TA, :DIC])
    ])

    # 5. Run the Round-Robin Iterations
    for (i, p1_sym) in enumerate(all_pars)
        # Look at the rest of the array after the current index
        for p2_sym in all_pars[i+1:end] 
            
            current_pair = Set([p1_sym, p2_sym])
            
            if current_pair in invalid_pairs
                continue
            end

            # Create a localized testset for every valid combination
            @testset "Input Pair: $p1_sym & $p2_sym" begin
                
                # Construct the keyword arguments dynamically
                val1 = calculated_vars[p1_sym]
                val2 = calculated_vars[p2_sym]
                input_kwargs = Dict(p1_sym => val1, p2_sym => val2)

                # Solve the system using the new pair
                rr_results = whole_system(; input_kwargs..., ENV_KWARGS...)

                # Assert that all major parameters match the baseline
                # Using an absolute tolerance of 1e-4
                tol = 1e-3

                @test isapprox(rr_results.TA, baseline.TA, atol=tol)
                @test isapprox(rr_results.DIC, baseline.DIC, atol=tol)
                @test isapprox(rr_results.pHtot, baseline.pHtot, atol=tol)
                @test isapprox(rr_results.pCO₂, baseline.pCO₂, atol=tol)
                @test isapprox(rr_results.CO₃, baseline.CO₃, atol=tol)
                @test isapprox(rr_results.HCO₃, baseline.HCO₃, atol=tol)
                @test isapprox(rr_results.CO₂, baseline.CO₂, atol=tol)
            end
        end
    end
end