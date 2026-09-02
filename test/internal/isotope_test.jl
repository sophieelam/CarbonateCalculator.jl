using Test
using CarbonateCalculator

# Individual isotope conversions, tested below the public API, so they come from the
# Isotopes submodule directly.
import CarbonateCalculator.Isotopes: get_alphaB, calc_ABT, ABOH3_from_H_ABT,
                                     ABOH4_from_H_ABT, A11_to_δ11, δ11_to_A11,
                                     get_ϵ, alphaB_to_ϵ, ϵ_to_alpha, calc_KB,
                                     pH_from_δ, pKB_from_δ, calc_δ11BT, calc_δ11B4, calc_ϵ

@testset "BoronIsotopeFunctions: Test B isotope functions" begin

    @testset "get_alphaB" begin
        @test get_alphaB() == boron_ref.alphaB
    end

    @testset "calc_ABT (H, BOH₃)" begin
        @test calc_ABT(H=boron_ref.H, ABOH₃=boron_ref.ABOH₃, Ks=boron_ref.Ks, alphaB=boron_ref.alphaB) ≈ boron_ref.ABT atol=1e-6
    end

    @testset "calc_ABT (H, BOH₄)" begin
        @test calc_ABT(H=boron_ref.H, ABOH₄=boron_ref.ABOH₄, Ks=boron_ref.Ks, alphaB=boron_ref.alphaB) ≈ boron_ref.ABT atol=1e-6
    end

    @testset "ABOH3_from_H_ABT" begin
        @test ABOH3_from_H_ABT(boron_ref.H, boron_ref.ABT, boron_ref.Ks, boron_ref.alphaB) ≈ boron_ref.ABOH₃ atol=1e-6
    end

    @testset "ABOH4_from_H_ABT" begin
        @test ABOH4_from_H_ABT(boron_ref.H, boron_ref.ABT, boron_ref.Ks, boron_ref.alphaB) ≈ boron_ref.ABOH₄ atol=1e-6
    end

    # Isotope unit conversions
    @testset "A11_to_δ11" begin
        @test A11_to_δ11(0.807817779214075) ≈ 39.5 atol=1e-6
    end

    @testset "δ11_to_A11" begin
        @test δ11_to_A11(39.5) ≈ 0.807817779214075 atol=1e-6
    end

end

# The δ-notation layer, which works in ‰ rather than fractional abundance. Tested by the
# identities that define it rather than against frozen numbers: each function has an inverse
# in the same set, so a round trip has to return the value it started from.
@testset "BoronIsotopeδ: the δ-notation layer" begin
    KB = boron_ref.Ks.KB
    pH = -log10(boron_ref.H)
    δBT = 39.61

    @testset "ϵ and alpha are inverses" begin
        @test ϵ_to_alpha(alphaB_to_ϵ(get_alphaB())) ≈ get_alphaB() atol=1e-12
        @test get_ϵ() ≈ (get_alphaB() - 1) * 1000 atol=1e-9
    end

    # calc_KB takes either species; supplying ABOH₄ must not discard it.
    #
    # Loose rtol because the reference abundances are stored to 8 decimal places and KB is
    # recovered from the *difference* between two values near 0.808, which is three orders
    # of magnitude smaller than they are. The rounding is what sets the accuracy here, not
    # the function: the exact round trips below hold to ~1e-11.
    @testset "calc_KB accepts either species" begin
        from_ABOH4 = calc_KB(boron_ref.H, boron_ref.alphaB, boron_ref.ABT, boron_ref.ABOH₄)
        from_ABOH3 = calc_KB(boron_ref.H, boron_ref.alphaB, boron_ref.ABT,
                             nothing, boron_ref.ABOH₃)
        @test from_ABOH4 ≈ boron_ref.Ks.KB rtol=1e-5
        @test from_ABOH4 ≈ from_ABOH3 rtol=1e-5
    end

    @testset "calc_δ11B4 and calc_δ11BT are inverses" begin
        δB4 = calc_δ11B4(pH, KB, δBT)
        @test calc_δ11BT(pH, KB, δB4) ≈ δBT atol=1e-8
    end

    @testset "pH_from_δ inverts calc_δ11B4" begin
        δB4 = calc_δ11B4(pH, KB, δBT)
        @test pH_from_δ(KB, δBT, δB4) ≈ pH atol=1e-8
    end

    @testset "pKB_from_δ recovers KB" begin
        δB4 = calc_δ11B4(pH, KB, δBT)
        @test pKB_from_δ(pH, δBT, δB4) ≈ -log10(KB) atol=1e-8
    end

    @testset "calc_ϵ recovers the fractionation used" begin
        δB4 = calc_δ11B4(pH, KB, δBT)
        @test calc_ϵ(pH, KB, δBT, δB4) ≈ get_ϵ() atol=1e-6
    end

    # An independent path to the same number: the public solver computes δBOH₄ through the
    # fractional-abundance route, so agreement checks the two against each other.
    @testset "agrees with the public solver" begin
        result = boron_isotopes(pHtot = 8.1, δBT = δBT, temp_c = 20.0, sal = 35.0)
        @test calc_δ11B4(8.1, result.Ks.KB, δBT) ≈ result.δBOH₄ rtol=1e-10
        @test pH_from_δ(result.Ks.KB, δBT, result.δBOH₄) ≈ 8.1 atol=1e-8
    end
end