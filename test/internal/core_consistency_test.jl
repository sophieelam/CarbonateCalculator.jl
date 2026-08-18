# The two cores must stay one implementation.
#
# `carbon_system_core` and `whole_system_core` were near-duplicates — 200 and 241 lines with
# 29 unique to the carbon one — and the duplication had produced two defects, each present
# in one copy and not the other. They now compose shared stages from `src/core.jl`; these
# tests are what stops them drifting apart again.

using Test
using CarbonateCalculator

@testset "Core consistency" begin

    @testset "the two systems agree on every shared output" begin
        # carbon_system's 40 output keys are a strict subset of whole_system's 50, so for a
        # boron-free input every shared key must match. This is the test that would catch a
        # stage being applied in one core and not the other.
        for kwargs in ((TA = 2300.0, DIC = 2000.0),
                       (TA = 2300.0, pHtot = 8.1),
                       (DIC = 2000.0, pCO₂ = 400.0),
                       (TA = 2300.0, DIC = 2000.0, temp_c = 2.0, pres_bar = 400.0),
                       (TA = 2300.0, DIC = 2000.0, PT = 2.0, SiT = 30.0,
                        H2ST = 5.0, NH4T = 3.0),
                       (TA = 2.3, DIC = 2.0, unit = "mmol"))
            carbon = carbon_system(; kwargs...)
            whole = whole_system(; kwargs...)
            for key in keys(carbon.val)
                value = carbon.val[key]
                value isa Real && isfinite(value) || continue
                @test whole.val[key] ≈ value rtol=1e-12
            end
        end
    end

    @testset "ΩA and ΩC are accepted as inputs by both" begin
        # whole_system used to throw `ArgumentError: Impossible! You haven't provided enough
        # information` here, because the Ω -> CO₃ input conversion existed only in the
        # carbon core. Ω is a statement about carbonate ion, so it is a valid second
        # parameter for either system.
        for (name, value) in ((:ΩC, 5.0), (:ΩA, 3.0))
            carbon = carbon_system(; TA = 2300.0, name => value)
            whole = whole_system(; TA = 2300.0, name => value)
            @test carbon.pHtot ≈ whole.pHtot rtol=1e-12
            # ...and it round-trips: what went in comes back out.
            @test getproperty(whole, name) ≈ value rtol=1e-10
            @test getproperty(carbon, name) ≈ value rtol=1e-10
        end
    end

    @testset "every concentration rescales with the unit" begin
        # `Alk_H2S` and `Alk_NH3` were missing from the carbon core's rescale list, so at
        # unit="mmol" they came back in mol/kg while every other concentration in the same
        # result was in mmol/kg — a factor of 1000 apart, inside one result.
        inputs = (TA = 2300.0, DIC = 2000.0, H2ST = 10.0, NH4T = 5.0)
        umol = carbon_system(; inputs..., unit = "umol")
        mmol = carbon_system(; TA = 2.3, DIC = 2.0, H2ST = 0.01, NH4T = 0.005,
                             unit = "mmol")

        # Same water, different units: every concentration should differ by exactly 1000.
        for key in (:Alk_H2S, :Alk_NH3, :CAlk, :BAlk, :OH, :HCO₃, :CO₃)
            @test umol.val[key] ≈ mmol.val[key] * 1e3 rtol=1e-9
        end

        # And the two systems must agree, which is what drifted before.
        for unit in ("umol", "mmol")
            scale = unit == "umol" ? 1.0 : 1e-3
            carbon = carbon_system(; TA = 2300.0 * scale, DIC = 2000.0 * scale,
                                   H2ST = 10.0 * scale, NH4T = 5.0 * scale, unit = unit)
            whole = whole_system(; TA = 2300.0 * scale, DIC = 2000.0 * scale,
                                 H2ST = 10.0 * scale, NH4T = 5.0 * scale, unit = unit)
            @test carbon.Alk_H2S ≈ whole.Alk_H2S rtol=1e-12
            @test carbon.Alk_NH3 ≈ whole.Alk_NH3 rtol=1e-12
        end
    end

end
