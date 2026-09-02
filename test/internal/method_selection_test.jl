# Parameterisation selection under `K_method = "automatic"`.
#
# `which_K` is no longer what the default reaches: the default is now KGen, which is a single
# parameterisation across the whole domain and so has none of the discontinuities that
# switching between fitted formulas introduces. `which_K` remains as the `"automatic"` choice,
# for callers who want the best-fitted formula for their conditions and accept the steps at the
# band edges that come with it.
#
# Within `which_K`, domain guards are tested before preference guards. Interleaving them —
# `temp_c > 35` ahead of the salinity guards — gives every sample above 35 °C Millero 2006
# whatever its salinity, and makes the freshwater and brine branches unreachable for exactly the
# warm estuaries and hypersaline lagoons they exist for.
#
# This file pins the selection table so that ordering cannot silently regress, and pins the
# split between the two paths.

using Test
using CarbonateCalculator
import CarbonateCalculator.Constants: which_K, MODERN_CALCIUM, MODERN_MAGNESIUM

@testset "K method selection" begin

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
        @test which_K(temp_c=25.0, sal=35.0) == "KGen"
    end

    @testset "ordinary seawater is unchanged" begin
        # The reordering must not touch the conditions almost every user is at. If any of
        # these move, the change was not confined to the domain guards.
        for temp_c in (0.0, 5.0, 15.0, 25.0, 35.0), sal in (30.0, 33.0, 35.0, 37.0)
            selected = which_K(; temp_c, sal)
            @test selected in ("KGen", "Millero 2002", "GP 1989")
        end
        @test which_K(temp_c=-1.0, sal=35.0) == "Millero 2002"
    end

    @testset "non-modern seawater still overrides everything" begin
        # Composition is the strongest signal: only KGen corrects for it, so it must win
        # over both the domain and the preference guards.
        @test which_K(temp_c=36.0, sal=55.0, Ca=0.02) == "KGen"
        @test which_K(temp_c=36.0, sal=0.5, Mg=0.03) == "KGen"
        # ...but an explicitly-modern composition is not "non-modern".
        @test which_K(temp_c=36.0, sal=0.5, Ca=MODERN_CALCIUM,
                      Mg=MODERN_MAGNESIUM) == "Millero 1979"
    end

    @testset "selection reaches the public API, and only through \"automatic\"" begin
        # A hot brine, where which_K picks Papadimitriou 2018 and so disagrees with KGen. That
        # disagreement is what makes this able to tell the two paths apart at all.
        conditions = (TA=2300.0, DIC=2000.0, temp_c=36.0, sal=55.0)

        automatic = carbon_system(; conditions..., K_method="automatic")
        selected = carbon_system(; conditions..., K_method="Papadimitriou 2018")
        @test automatic.pHtot ≈ selected.pHtot atol=1e-12

        # The default does *not* consult which_K any more — it is KGen everywhere, which is
        # the point: one parameterisation across the domain has no band edges to step across.
        default = carbon_system(; conditions...)
        @test default.pHtot ≈ carbon_system(; conditions..., K_method="KGen").pHtot atol=1e-12
        @test !isapprox(default.pHtot, automatic.pHtot; atol=1e-6)
    end

    @testset "the default is continuous where which_K is not" begin
        # which_K switches parameterisation at 2 °C for ordinary salinity, which moves pH by a
        # step far larger than the interval producing it. Anything differentiating with respect
        # to temperature across that point gets a meaningless answer, which is what choosing a
        # single parameterisation as the default avoids.
        pH(t; kw...) = carbon_system(; TA=2300.0, DIC=2000.0, sal=35.0, temp_c=t, kw...).pHtot

        step_default = abs(pH(2.0) - pH(1.99))
        step_automatic = abs(pH(2.0; K_method="automatic") - pH(1.99; K_method="automatic"))

        @test step_default < 1e-3           # ~1.7e-4, the ordinary 0.01 °C change
        @test step_automatic > 1e-2         # ~1.7e-2, a hundredfold jump at the band edge
    end

end
