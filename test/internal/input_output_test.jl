const δBT_SW = 39.61  # modern seawater δ¹¹B, matches Isotopes.get_δBT()

@testset "TestInputOutput: Test internal consistency of measured/target condition calculations" begin

    # A sample measured at one temperature/pressure and reported at another must survive the
    # round trip unchanged: recalculate to the collection conditions, treat that state as a
    # fresh measurement, and come back. That is the route a caller takes when re-deriving bench
    # values from reported ones.
    #
    # There was a second, shorter route here — carry the result straight back — which is no
    # longer possible: a collection state cannot be carried again, deliberately, because a
    # sample has one set of collection conditions. Going through a declared fresh measurement
    # is what replaces it, and tests the same inversion with one more step.

    @testset "Csys" begin

        @testset "Temperature Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, temp_c=20)
            at_target = at_collection_conditions(measured, temp_c=30)

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, temp_c=30)
            @test measured.pHtot ≈ at_collection_conditions(remeasured, temp_c=20).pHtot atol=1e-6
        end

        @testset "Pressure Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, pres_bar=0)
            at_target = at_collection_conditions(measured, pres_bar=400)

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, pres_bar=400)
            @test measured.pHtot ≈ at_collection_conditions(remeasured, pres_bar=0).pHtot atol=1e-6
        end

        @testset "Temperature and Pressure together" begin
            measured = whole_system(pHtot=8.1, TA=2300, temp_c=20, pres_bar=0)
            at_target = at_collection_conditions(measured, temp_c=2, pres_bar=400)

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, temp_c=2, pres_bar=400)
            @test measured.pHtot ≈ at_collection_conditions(remeasured, temp_c=20,
                                                            pres_bar=0).pHtot atol=1e-6
        end

    end

    @testset "CBsys" begin

        @testset "Temperature Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, temp_c=20, δBT=δBT_SW)
            at_target = at_collection_conditions(measured, temp_c=30)

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, temp_c=30, δBT=δBT_SW)
            @test measured.pHtot ≈ at_collection_conditions(remeasured, temp_c=20).pHtot atol=1e-6
        end

        @testset "Pressure Effect" begin
            measured = whole_system(pHtot=8.1, TA=2300, pres_bar=0, δBT=δBT_SW)
            at_target = at_collection_conditions(measured, pres_bar=400)

            remeasured = whole_system(pHtot=at_target.pHtot, TA=2300, pres_bar=400, δBT=δBT_SW)
            @test measured.pHtot ≈ at_collection_conditions(remeasured, pres_bar=0).pHtot atol=1e-6
        end

    end

    @testset "Retired output-condition arguments" begin
        # These names selected a second set of conditions in the two-condition API. Because
        # every calculation function ends in `kwargs...`, an unrecognised name is otherwise
        # absorbed silently and the call computes at the default 25 °C.
        @test_throws ArgumentError whole_system(pHtot=8.1, TA=2300, temp_c=20, T_out=30)
        @test_throws ArgumentError whole_system(pHtot=8.1, TA=2300, pres_bar=0, P_out=400)

        # Salinity was never a physical degree of freedom here: a water sample's salinity
        # does not change between collection and measurement.
        @test_throws ArgumentError whole_system(pHtot=8.1, TA=2300, sal=28.2, S_out=38.1)
    end

end
