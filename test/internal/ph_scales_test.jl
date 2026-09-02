# The carbon and boron entry points convert between pH scales through different code:
# carbon_system_core derives them inline, while boron_system and boron_isotopes go through
# Helpers.calc_pH_scale. Nothing checked that the two agreed, and they did not — calc_pH_scale
# had the sign of the SWS-to-total term wrong in its NBS offset, putting pHNBS out by twice
# that offset (0.02 pH at S=35) on the boron path only.

@testset "pH scale conversions" begin

    @testset "carbon and boron paths agree" begin
        for temp_c in (2.0, 15.0, 25.0, 35.0), sal in (30.0, 35.0, 38.0)
            b = boron_system(; pHtot=8.1, BT=415.7, δBT=39.61, temp_c=temp_c, sal=sal)
            c = carbon_system(; pHtot=8.1, TA=2300.0, temp_c=temp_c, sal=sal).val

            @test b.pHtot  ≈ c.pHtot  atol=1e-12
            @test b.pHsws  ≈ c.pHsws  atol=1e-12
            @test b.pHfree ≈ c.pHfree atol=1e-12
            @test b.pHNBS  ≈ c.pHNBS  atol=1e-12
        end
    end

    @testset "NBS offset follows the definition of fH" begin
        # fH = a(H)_NBS / m(H)_SWS, so pHNBS = pHsws - log10(fH). Checking against the
        # definition rather than against the other implementation, so that a matching pair
        # of wrong implementations cannot pass.
        for temp_c in (2.0, 25.0), sal in (30.0, 35.0)
            c = carbon_system(; pHtot=8.1, TA=2300.0, temp_c=temp_c, sal=sal).val
            fH = CarbonateCalculator.Helpers.calc_fH(temp_c, sal)
            @test c.pHNBS ≈ c.pHsws - log10(fH) atol=1e-10
        end
    end

    @testset "fH is physically plausible" begin
        # ~0.71 for seawater. Feeding Kelvin to the polynomial instead of Celsius returns
        # ~1.76, which this range rejects.
        for temp_c in (0.0, 25.0, 40.0), sal in (30.0, 35.0, 40.0)
            @test 0.6 < CarbonateCalculator.Helpers.calc_fH(temp_c, sal) < 0.85
        end
        @test CarbonateCalculator.Helpers.calc_fH(25.0, 35.0) ≈ 0.71340431875 rtol=1e-12
    end

end
