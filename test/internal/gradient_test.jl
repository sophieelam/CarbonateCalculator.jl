# Buffer factors computed on demand: ∂X/∂Y at constant Z.
#
# The two hand-written factors in carbon.jl are the reference. They are called by the pipeline
# *before* the concentrations are rescaled, so they work in mol/kg — which is why the reference
# result below is built with `unit = "mol"` rather than passing values off a µmol/kg result.
# Doing the latter returns a plausible number about 1% wrong.

using Test
using CarbonateCalculator

@testset "Gradients and buffer factors" begin

    kw = (TA = 2300.0, DIC = 2000.0, temp_c = 20.0, sal = 35.0, PT = 0.5, SiT = 3.0)
    result = carbon_system(; kw...)

    @testset "agrees with the hand-written factors" begin
        in_mol = carbon_system(; TA = 2300e-6, DIC = 2000e-6, temp_c = 20.0, sal = 35.0,
                               PT = 0.5e-6, SiT = 3e-6, unit = "mol")
        v = in_mol.val

        reference = CarbonateCalculator.calc_revelle_factor(
            v.TA, v.DIC, v.BT, v.PT, v.SiT, v.ST, v.FT, v.H2ST, v.NH4T, in_mol.Ks)
        @test revelle_factor(in_mol) ≈ reference rtol=1e-12

        capacity = CarbonateCalculator.calc_buffer_capacity(
            v.pHtot, v.DIC, v.BT, v.PT, v.SiT, v.ST, v.FT, v.H2ST, v.NH4T, in_mol.Ks)
        @test calc_gradient(in_mol, :TA, :pHtot; constant = :DIC) ≈ capacity rtol=1e-12
    end

    @testset "the held constant is half the definition" begin
        at_TA = calc_gradient(result, :pCO₂, :DIC; constant = :TA)
        at_pH = calc_gradient(result, :pCO₂, :DIC; constant = :pHtot)

        # Not a small difference: roughly a factor of ten. Naming numerator and denominator
        # alone would not have specified which of these was meant.
        @test at_TA > 5 * at_pH
        @test at_TA > 0 && at_pH > 0
    end

    @testset "relative is the gradient normalised" begin
        gradient = calc_gradient(result, :pCO₂, :DIC; constant = :TA)
        relative = calc_relative_gradient(result, :pCO₂, :DIC; constant = :TA)
        @test relative ≈ gradient * result.DIC / result.pCO₂ rtol=1e-14

        # fCO₂ and pCO₂ differ by a factor fixed at constant T, S and P, so it cancels.
        @test revelle_factor(result) ≈
              calc_relative_gradient(result, :pCO₂, :DIC; constant = :TA) rtol=1e-14

        # Dimensionless, so the reporting unit cannot matter.
        in_mol = carbon_system(; TA = 2300e-6, DIC = 2000e-6, temp_c = 20.0, sal = 35.0,
                               PT = 0.5e-6, SiT = 3e-6, unit = "mol")
        @test revelle_factor(result) ≈ revelle_factor(in_mol) rtol=1e-10
    end

    @testset "ill-posed requests are refused" begin
        # Varying and holding the same quantity.
        @test_throws "both varied and held" calc_gradient(result, :pCO₂, :DIC; constant = :DIC)

        # pCO₂ and fCO₂ are one quantity expressed two ways, so they cannot be an independent
        # pair. This comes from the solver's own determinacy check rather than a second list.
        @test_throws ArgumentError calc_gradient(result, :pHtot, :fCO₂; constant = :pCO₂)

        # A name that is not a computed value at all.
        @test_throws Exception calc_gradient(result, :pCO₂, :DIC; constant = :not_a_parameter)
    end

    @testset "works for every kind of result" begin
        # Boron in scope: the constraint set comes from PARAMETER_GROUPS, so the boron
        # parameters are stripped too and cannot over-determine the re-solve.
        whole = whole_system(; kw..., δBT = 39.61)
        @test revelle_factor(whole) > 0
        @test calc_gradient(whole, :pHtot, :DIC; constant = :TA) < 0   # more DIC, lower pH

        # A result already carried to its collection conditions.
        collected = at_collection_conditions(result; temp_c = 2.0)
        @test revelle_factor(collected) != revelle_factor(result)

        # A result whose measured pair was not TA/DIC, so the held value is one it derived.
        from_pH = carbon_system(; pHtot = 8.1, TA = 2300.0, temp_c = 20.0, sal = 35.0)
        @test revelle_factor(from_pH) > 0
    end

    @testset "with_gradient adds a value without disturbing the result" begin
        added = with_gradient(result, :revelle, :fCO₂, :DIC; constant = :TA, relative = true)

        @test added.revelle ≈ revelle_factor(result) rtol=1e-14
        @test keys(added.val) == (keys(result.val)..., :revelle)
        @test !haskey(result.val, :revelle)          # the original is untouched
        @test added.pHtot == result.pHtot            # everything else survives
        @test added.system == result.system
        @test added.is_collection_state == result.is_collection_state

        # Absolute as well as relative.
        raw = with_gradient(result, :dpCO2_dDIC, :pCO₂, :DIC; constant = :TA)
        @test raw.dpCO2_dDIC ≈ calc_gradient(result, :pCO₂, :DIC; constant = :TA) rtol=1e-14

        # A name already in use would silently shadow a computed value.
        @test_throws "already has a value called" with_gradient(
            result, :pHtot, :pCO₂, :DIC; constant = :TA)

        # The flag survives, so a carried result keeps refusing a second carry.
        collected = with_gradient(at_collection_conditions(result; temp_c = 2.0),
                                  :revelle, :fCO₂, :DIC; constant = :TA, relative = true)
        @test collected.is_collection_state
        @test_throws ArgumentError at_collection_conditions(collected; temp_c = 4.0)
    end

end
