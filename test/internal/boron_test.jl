using Test
using CarbonateCalculator

# Individual speciation equations, tested below the public API, so they come from the
# Boron submodule directly.
import CarbonateCalculator.Boron: H_from_BT_BOH3, H_from_BT_BOH4, BT_from_pH_BOH3,
                                  BT_from_pH_BOH4, calc_BOH3, calc_BOH4, calc_chiB

@testset "BoronFunctions: Test B concentration and speciation functions" begin

    @testset "H_from_BT_BOH3" begin
        @test H_from_BT_BOH3(boron_ref.BT, boron_ref.BOH₃, boron_ref.Ks) ≈ boron_ref.H atol=1e-6
    end

    @testset "H_from_BT_BOH4" begin
        @test H_from_BT_BOH4(boron_ref.BT, boron_ref.BOH₄, boron_ref.Ks) ≈ boron_ref.H atol=1e-6
    end

    @testset "BT_from_pH_BOH3" begin
        @test BT_from_pH_BOH3(boron_ref.pHtot, boron_ref.BOH₃, boron_ref.Ks) ≈ boron_ref.BT atol=1e-6
    end

    @testset "BT_from_pH_BOH4" begin
        @test BT_from_pH_BOH4(boron_ref.pHtot, boron_ref.BOH₄, boron_ref.Ks) ≈ boron_ref.BT atol=1e-6
    end

    @testset "calc_BOH3" begin
        @test calc_BOH3(boron_ref.BT, boron_ref.H, boron_ref.Ks) ≈ boron_ref.BOH₃ atol=1e-6
    end

    @testset "calc_BOH4" begin
        @test calc_BOH4(boron_ref.BT, boron_ref.H, boron_ref.Ks) ≈ boron_ref.BOH₄ atol=1e-6
    end

    @testset "calc_chiB" begin
        @test calc_chiB(boron_ref.H, boron_ref.Ks) == 1 / (1 + boron_ref.Ks.KB / boron_ref.H)
    end

end