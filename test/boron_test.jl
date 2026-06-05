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