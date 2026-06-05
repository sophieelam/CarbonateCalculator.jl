using .Calculator
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