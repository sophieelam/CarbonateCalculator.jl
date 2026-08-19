# Default parameterisation selection.
#
# `which_K` tests domain guards before preference guards. Interleaving them — `temp_c > 35`
# ahead of the salinity guards — gives every sample above 35 °C Millero 2006 whatever its
# salinity, and makes the freshwater and brine branches unreachable for exactly the warm
# estuaries and hypersaline lagoons they exist for.
#
# This file pins the selection table so that ordering cannot silently regress.

using Test
using CarbonateCalculator
import CarbonateCalculator.Constants: which_K, MODERN_CALCIUM, MODERN_MAGNESIUM

@testset "Default K method selection" begin

    @testset "salinity domain beats temperature preference" begin
        # The regression cases. Each was Millero 2006 before, outside the salinity range it
        # was fitted over.
        @test which_K(temp_c=36.0, sal=0.5) == "Millero 1979"
        @test which_K(temp_c=36.0, sal=55.0) == "Papadimitriou 2018"
        @test which_K(temp_c=45.0, sal=0.0) == "Millero 1979"
        @test which_K(temp_c=45.0, sal=60.0) == "Papadimitriou 2018"

        # ...and the same domains are still honoured below 35 °C, where they always were.
        @test which_K(temp_c=10.0, sal=0.5) == "Millero 1979"
        @test which_K(temp_c=10.0, sal=55.0) == "Papadimitriou 2018"
    end

    @testset "temperature preference within the normal salinity range" begin
        @test which_K(temp_c=36.0, sal=35.0) == "Millero 2006"
        @test which_K(temp_c=1.0, sal=35.0) == "Millero 2002"
        @test which_K(temp_c=-1.0, sal=20.0) == "GP 1989"
        @test which_K(temp_c=25.0, sal=35.0) == "MyAMI"
    end

    @testset "ordinary seawater is unchanged" begin
        # The reordering must not touch the conditions almost every user is at. If any of
        # these move, the change was not confined to the domain guards.
        for temp_c in (0.0, 5.0, 15.0, 25.0, 35.0), sal in (30.0, 33.0, 35.0, 37.0)
            selected = which_K(; temp_c, sal)
            @test selected in ("MyAMI", "Millero 2002", "GP 1989")
        end
        @test which_K(temp_c=-1.0, sal=35.0) == "Millero 2002"
    end

    @testset "non-modern seawater still overrides everything" begin
        # Composition is the strongest signal: only MyAMI corrects for it, so it must win
        # over both the domain and the preference guards.
        @test which_K(temp_c=36.0, sal=55.0, Ca=0.02) == "MyAMI"
        @test which_K(temp_c=36.0, sal=0.5, Mg=0.03) == "MyAMI"
        # ...but an explicitly-modern composition is not "non-modern".
        @test which_K(temp_c=36.0, sal=0.5, Ca=MODERN_CALCIUM,
                      Mg=MODERN_MAGNESIUM) == "Millero 1979"
    end

    @testset "selection reaches the public API" begin
        # which_K is only useful if K_method="default" actually consults it, so check the
        # selected method and the default agree at a condition where the fix applies.
        hot_brine = carbon_system(TA=2300.0, DIC=2000.0, temp_c=36.0, sal=55.0)
        pinned = carbon_system(TA=2300.0, DIC=2000.0, temp_c=36.0, sal=55.0,
                               K_method="Papadimitriou 2018")
        @test hot_brine.pHtot ≈ pinned.pHtot atol=1e-12
    end

end
