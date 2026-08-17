const δBT_SW = 39.61  # modern seawater δ¹¹B, matches Isotopes.get_δBT()

@testset "TestInputOutput: Test internal consistency of measured/target condition calculations" begin

    # A sample measured at one temperature/pressure and reported at another must survive the
    # round trip unchanged. Two routes are checked because they can fail independently:
    #
    #   direct   - recalculate to the target, then recalculate straight back
    #   remeasure - recalculate to the target, treat that as a fresh measurement, come back
    #
    # The second is the one the old two-condition API exercised.

    @testset "Csys" begin

        @testset "Temperature Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, temp_c=20)
            at_target = recalculate_at_target_conditions(measured, temp_c=30)

            back = recalculate_at_target_conditions(at_target, temp_c=20)
            @test measured.pHtot ≈ back.pHtot atol=1e-6

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, temp_c=30)
            @test measured.pHtot ≈ recalculate_at_target_conditions(remeasured, temp_c=20).pHtot atol=1e-6
        end

        @testset "Pressure Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, pres_bar=0)
            at_target = recalculate_at_target_conditions(measured, pres_bar=400)

            back = recalculate_at_target_conditions(at_target, pres_bar=0)
            @test measured.pHtot ≈ back.pHtot atol=1e-6

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, pres_bar=400)
            @test measured.pHtot ≈ recalculate_at_target_conditions(remeasured, pres_bar=0).pHtot atol=1e-6
        end

        @testset "Temperature and Pressure together" begin
            measured = whole_system(pHtot=8.1, TA=2300, temp_c=20, pres_bar=0)
            at_target = recalculate_at_target_conditions(measured, temp_c=2, pres_bar=400)
            back = recalculate_at_target_conditions(at_target, temp_c=20, pres_bar=0)
            @test measured.pHtot ≈ back.pHtot atol=1e-6
        end

    end

    @testset "CBsys" begin

        @testset "Temperature Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, temp_c=20, δBT=δBT_SW)
            at_target = recalculate_at_target_conditions(measured, temp_c=30)

            back = recalculate_at_target_conditions(at_target, temp_c=20)
            @test measured.pHtot ≈ back.pHtot atol=1e-6

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, temp_c=30, δBT=δBT_SW)
            @test measured.pHtot ≈ recalculate_at_target_conditions(remeasured, temp_c=20).pHtot atol=1e-6
        end

        @testset "Pressure Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, pres_bar=0, δBT=δBT_SW)
            at_target = recalculate_at_target_conditions(measured, pres_bar=400)

            back = recalculate_at_target_conditions(at_target, pres_bar=0)
            @test measured.pHtot ≈ back.pHtot atol=1e-6

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, pres_bar=400, δBT=δBT_SW)
            @test measured.pHtot ≈ recalculate_at_target_conditions(remeasured, pres_bar=0).pHtot atol=1e-6
        end

    end

    @testset "Retired output-condition arguments" begin
        # These used to select a second set of conditions. They are gone, and because every
        # calculation function ends in `kwargs...` an unrecognised name would otherwise be
        # absorbed silently and the call would compute at the default 25 °C.
        @test_throws ArgumentError whole_system(pHtot=8.1, TA=2300, temp_c=20, T_out=30)
        @test_throws ArgumentError whole_system(pHtot=8.1, TA=2300, pres_bar=0, P_out=400)

        # Salinity was never a physical degree of freedom here: a water sample's salinity
        # does not change between collection and measurement.
        @test_throws ArgumentError whole_system(pHtot=8.1, TA=2300, sal=28.2, S_out=38.1)
    end

end
