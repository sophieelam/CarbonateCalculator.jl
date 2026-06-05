using Random
using Test # Assuming you have this imported elsewhere, but good to ensure!

Random.seed!(42)

n = 101
# Generate the arrays (adding _arr to the names makes it clear they are vectors)
pHtot_arr = 7.6 .+ rand(n) .* (8.8 - 7.6)
DIC_arr = 1900.0 .+ rand(n) .* (2200.0 - 1900.0)
BT_arr = 350.0 .+ rand(n) .* (450.0 - 350.0)
δBT_arr = 35.0 .+ rand(n) .* (45.0 - 35.0)
T_arr = 15.0 .+ rand(n) .* (35.0 - 15.0)
S_arr = 30.0 .+ rand(n) .* (40.0 - 30.0)

@testset "ReferenceDataTestCase: Test boron isotopes" begin

    @testset "Bisotopes" begin
        
        # Loop through all 101 random values and test them as scalars
        for i in 1:n
            pHtot = pHtot_arr[i]
            DIC = DIC_arr[i]
            BT = BT_arr[i]
            δBT = δBT_arr[i]
            T = T_arr[i]
            S = S_arr[i]

            # Generate the base test system for THIS specific iteration
            test_sys = whole_system(pHtot=pHtot, δBT=δBT, BT=BT, DIC=DIC, T_in=T, S_in=S)

            # Test boron_isotopes
            check_A = boron_isotopes(δBOH₄=test_sys.δBOH₄, δBT=test_sys.δBT, T_in=T, S_in=S)
            @test test_sys.pHtot ≈ check_A.pHtot rtol=1e-10
            
            # Test boron_system
            check_B = boron_system(δBOH₄=test_sys.δBOH₄, δBT=test_sys.δBT, BT=test_sys.BT, T_in=T, S_in=S)
            @test test_sys.pHtot ≈ check_B.pHtot rtol=1e-10

            # Test whole_system
            check_CB = whole_system(δBOH₄=test_sys.δBOH₄, δBT=test_sys.δBT, DIC=test_sys.DIC, BT=test_sys.BT, T_in=T, S_in=S)
            @test test_sys.pHtot ≈ check_CB.pHtot rtol=1e-10
        end
        
    end

end
