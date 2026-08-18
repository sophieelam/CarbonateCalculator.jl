using Test
using CarbonateCalculator

# Individual speciation equations, tested below the public API, so they come from the
# Carbon submodule directly.
import CarbonateCalculator.Carbon: CO₂_from_pH_DIC, H_from_HCO₃_CO₃, H_from_HCO₃_TA,
                                   pH_from_HCO₃_DIC, H_from_CO₃_TA, H_from_CO₃_DIC,
                                   pH_from_TA_DIC, calc_CO₂, calc_CO₃, calc_HCO₃, calc_TA, calc_TA_components,
                                   fCO₂_to_CO₂, CO₂_to_fCO₂, fCO₂_to_pCO₂, pCO₂_to_fCO₂,
                                   DIC_from_CO₂_pH, H_from_CO₂_HCO₃, H_from_CO₂_CO₃,
                                   pH_from_CO₂_TA, H_from_CO₂_DIC, DIC_from_pH_HCO₃,
                                   DIC_from_pH_CO₃, DIC_from_pH_TA

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
        TA, CAlk, BAlk, PAlk, SiAlk, OH, Hfree, HSO₄, HF, Alk_H2S, Alk_NH3 = calc_TA_components(
            carbon_ref.H,
            carbon_ref.DIC / carbon_ref.unit,
            carbon_ref.BT / carbon_ref.unit,
            carbon_ref.TP / carbon_ref.unit,
            carbon_ref.TSi / carbon_ref.unit,
            carbon_ref.TS,
            carbon_ref.TF,
            0.0, # H2ST
            0.0, # NH4T
            carbon_ref.Ks
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


@testset "iterative solves converge away from pH 8" begin
    # Every iterative solve used to start at pH 8.0 and stay there, so an answer far from
    # seawater diverged: TA+DIC raised ConvergenceFailed above pH 10.5, CO₂+TA above pH 9.5,
    # HCO₃+TA above pH 10. Nothing caught it, because every existing case sits near pH 8.
    #
    # Each case here is generated from (pH, DIC), which is analytic and unambiguous, and then
    # solved back from the derived pair. The three pairs below have residuals that are
    # monotonic in pH, so there is exactly one root to find and the expected answer is not in
    # question. CO₃+TA and HCO₃+DIC are deliberately absent: both are genuinely ambiguous over
    # part of this range and are recorded as open items rather than pinned here.
    env = CarbonateCalculator.Constants.calculate_constants(temp_c = 25.0, sal = 35.0)
    Ks = env.Ks
    totals = (BT = env.BT, PT = 0.0, SiT = 0.0, ST = env.ST, FT = env.FT,
              H2ST = 0.0, NH4T = 0.0)

    "Everything implied by a (pH, DIC) pair, in mol/kg."
    function reference(pH, DIC)
        H = 10.0^(-pH)
        TA = calc_TA(H, DIC, totals.BT, totals.PT, totals.SiT, totals.ST, totals.FT,
                      totals.H2ST, totals.NH4T, Ks)
        return (TA = TA, CO₂ = calc_CO₂(H, DIC, Ks), HCO₃ = calc_HCO₃(H, DIC, Ks))
    end

    # DIC low enough that the corresponding pH is far above the old fixed starting point.
    for (pH, DIC) in ((9.5, 500e-6), (10.0, 500e-6), (10.5, 2000e-6), (11.0, 200e-6))
        state = reference(pH, DIC)

        @test pH_from_CO₂_TA(state.CO₂, state.TA, totals.BT, totals.PT, totals.SiT,
                             totals.ST, totals.FT, totals.H2ST, totals.NH4T, Ks) ≈ pH atol=1e-6

        @test pH_from_TA_DIC(state.TA, DIC, totals.BT, totals.PT, totals.SiT, totals.ST,
                             totals.FT, totals.H2ST, totals.NH4T, Ks) ≈ pH atol=1e-6

        H = H_from_HCO₃_TA(state.HCO₃, state.TA, totals.BT, totals.PT, totals.SiT, totals.ST,
                           totals.FT, totals.H2ST, totals.NH4T, Ks)
        @test -log10(H) ≈ pH atol=1e-3
    end

    # The fallback must not cost the partials. A bracketing `find_zero` would return a plain
    # Float64 here and σ would come back as exactly zero — see `_solve_pH`.
    solve = CarbonateSystem(:carbon; varying = (:TA, :DIC), varying_errors = (:TA, :DIC))
    high_pH = solve(2300.0, 180.0, 2.0, 2.0)
    @test high_pH.pHtot > 10.0
    @test high_pH.err.pHtot > 0.0
end