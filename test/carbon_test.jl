@testset "CarbonFnTestCase: Test all C functions" begin

    @testset "DIC_from_CO₂_pH" begin
        @test DIC_from_CO₂_pH(carbon_ref.CO₂, carbon_ref.pHtot, carbon_ref.Ks) ≈ carbon_ref.DIC atol=1e-6
    end

    @testset "H_from_CO₂_HCO₃ (zf)" begin
        @test H_from_CO₂_HCO₃(carbon_ref.CO₂, carbon_ref.HCO₃, carbon_ref.Ks) ≈ carbon_ref.H atol=1e-6
    end

    @testset "H_from_CO₂_CO₃ (zf)" begin
        @test H_from_CO₂_CO₃(carbon_ref.CO₂, carbon_ref.CO₃, carbon_ref.Ks) ≈ carbon_ref.H atol=1e-6
    end

    @testset "pH_from_CO₂_TA" begin
        res = pH_from_CO₂_TA(
            carbon_ref.CO₂ / carbon_ref.unit,
            carbon_ref.TA / carbon_ref.unit,
            carbon_ref.BT / carbon_ref.unit,
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks
        )[1]
        @test res ≈ carbon_ref.pHtot atol=1e-6
    end

    @testset "H_from_CO₂_DIC (zf)" begin
        @test H_from_CO₂_DIC(carbon_ref.CO₂, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.H atol=1e-6
    end

    @testset "DIC_from_pH_HCO₃" begin
        @test DIC_from_pH_HCO₃(carbon_ref.pHtot, carbon_ref.HCO₃, carbon_ref.Ks) ≈ carbon_ref.DIC atol=1e-6
    end

    @testset "DIC_from_pH_CO₃" begin
        @test DIC_from_pH_CO₃(carbon_ref.pHtot, carbon_ref.CO₃, carbon_ref.Ks) ≈ carbon_ref.DIC atol=1e-6
    end

    @testset "DIC_from_pH_TA" begin
        res = DIC_from_pH_TA(
            carbon_ref.pHtot,
            carbon_ref.TA / carbon_ref.unit,
            carbon_ref.BT / carbon_ref.unit,
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks
        ) * carbon_ref.unit
        @test res ≈ carbon_ref.DIC atol=1e-6
    end

    @testset "CO₂_from_pH_DIC" begin
        @test CO₂_from_pH_DIC(carbon_ref.pHtot, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.CO₂ atol=1e-6
    end

    @testset "H_from_HCO₃_CO₃ (zf)" begin
        @test H_from_HCO₃_CO₃(carbon_ref.HCO₃, carbon_ref.CO₃, carbon_ref.Ks) ≈ carbon_ref.H atol=1e-6
    end

    @testset "H_from_HCO₃_TA (zf)" begin
        @test H_from_HCO₃_TA(
            carbon_ref.HCO₃ / carbon_ref.unit, 
            carbon_ref.TA / carbon_ref.unit, 
            carbon_ref.BT / carbon_ref.unit, 
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks
        ) ≈ carbon_ref.H atol=1e-6
    end

    @testset "pH_from_HCO₃_DIC (zf)" begin
        @test pH_from_HCO₃_DIC(carbon_ref.HCO₃, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.pHtot atol=1e-6
    end

    @testset "H_from_CO₃_TA (zf)" begin
        @test H_from_CO₃_TA(
            carbon_ref.CO₃ / carbon_ref.unit, 
            carbon_ref.TA / carbon_ref.unit, 
            carbon_ref.BT / carbon_ref.unit, 
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks
        ) ≈ carbon_ref.H atol=1e-6
    end

    @testset "H_from_CO₃_DIC (zf)" begin
        @test H_from_CO₃_DIC(carbon_ref.CO₃, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.H atol=1e-6
    end

    @testset "pH_from_TA_DIC" begin
        res = pH_from_TA_DIC(
            carbon_ref.TA / carbon_ref.unit,
            carbon_ref.DIC / carbon_ref.unit,
            carbon_ref.BT / carbon_ref.unit,
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks
        )[1]
        @test res ≈ carbon_ref.pHtot atol=1e-6
    end

    @testset "calc_CO₂" begin
        @test calc_CO₂(carbon_ref.H, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.CO₂ atol=1e-6
    end

    @testset "calc_CO₃" begin
        @test calc_CO₃(carbon_ref.H, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.CO₃ atol=1e-6
    end

    @testset "calc_HCO₃" begin
        @test calc_HCO₃(carbon_ref.H, carbon_ref.DIC, carbon_ref.Ks) ≈ carbon_ref.HCO₃ atol=1e-6
    end

    @testset "calc_TA" begin
        TA, CAlk, BAlk, PAlk, SiAlk, OH, Hfree, HSO₄, HF, Alk_H2S, Alk_NH3 = calc_TA(
            carbon_ref.H,
            carbon_ref.DIC / carbon_ref.unit,
            carbon_ref.BT / carbon_ref.unit,
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks;
            mode="multi" # Pass this as a keyword argument now!
        )

        @test TA * carbon_ref.unit ≈ carbon_ref.TA atol=1e-6
        @test CAlk * carbon_ref.unit ≈ carbon_ref.CAlk atol=1e-6
        @test BAlk * carbon_ref.unit ≈ carbon_ref.BAlk atol=1e-6
        @test PAlk * carbon_ref.unit ≈ carbon_ref.PAlk atol=1e-6
        @test SiAlk * carbon_ref.unit ≈ carbon_ref.SiAlk atol=1e-6
        @test OH * carbon_ref.unit ≈ carbon_ref.OH atol=1e-6
        @test Hfree * carbon_ref.unit ≈ carbon_ref.Hfree atol=1e-6
        @test HSO₄ * carbon_ref.unit ≈ carbon_ref.HSO₄ atol=1e-6
        @test HF * carbon_ref.unit ≈ carbon_ref.HF atol=1e-6
    end

    @testset "fCO₂_to_CO₂" begin
        @test fCO₂_to_CO₂(carbon_ref.fCO₂, carbon_ref.Ks) ≈ carbon_ref.CO₂ atol=1e-6
    end

    @testset "CO₂_to_fCO₂" begin
        @test CO₂_to_fCO₂(carbon_ref.CO₂, carbon_ref.Ks) ≈ carbon_ref.fCO₂ atol=1e-6
    end

    @testset "fCO₂_to_pCO₂" begin
        @test fCO₂_to_pCO₂(carbon_ref.fCO₂, carbon_ref.T) ≈ carbon_ref.pCO₂ atol=1e-5
    end

    @testset "pCO₂_to_fCO₂" begin
        @test pCO₂_to_fCO₂(carbon_ref.pCO₂, carbon_ref.T) ≈ carbon_ref.fCO₂ atol=1e-5
    end

end