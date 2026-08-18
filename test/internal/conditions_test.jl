# Tests for the two-stage conditions API: carbon_system/whole_system solve at the conditions
# a sample was measured at, recalculate_at_target_conditions re-solves at the conditions it
# was collected at.

const δBT_MODERN = 39.61

@testset "Two-stage conditions API" begin

    measured_kw = (TA=2300.0, DIC=2000.0, temp_c=20.0, sal=35.0, pres_bar=0.0)

    @testset "No-op recalculation reproduces the measured state" begin
        # The most important test here. Stage 2 rebuilds its call from what the result
        # carries, so anything it fails to carry shows up as a state that does not
        # reproduce. Run across the method matrix, not just defaults - the whole class of
        # bug this refactor fixes was a setting being dropped between conditions.
        k_methods = ["Roy 1993", "GP 1989", "Hansson 1973", "DM 1987", "HM 1973",
                     "Mehrbach 1973 A", "Mehrbach 1973 B", "CW 2003", "Lueker 2000",
                     "MPM 2002", "Millero 2002", "Millero 2006", "Millero 2010",
                     "Waters 2014", "SB 2020", "Sulpis 2020", "MyAMI"]

        for k in k_methods
            m = carbon_system(; measured_kw..., K_method=k)
            r = recalculate_at_target_conditions(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.CO₃ ≈ m.val.CO₃ rtol=1e-12
            @test r.val.ΩC ≈ m.val.ΩC rtol=1e-12
        end

        for kso4 in ["default", "Khoo", "WM13"], kf in ["default", "Perez"]
            m = carbon_system(; measured_kw..., K_method="Lueker 2000",
                              KSO4_method=kso4, KF_method=kf)
            r = recalculate_at_target_conditions(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.pHfree ≈ m.val.pHfree rtol=1e-12
        end

        for bt in ["default", "Lee", "KSK18"], ca in ["default", "Culkin", "RT67"]
            m = carbon_system(; measured_kw..., BT_method=bt, Ca_method=ca)
            r = recalculate_at_target_conditions(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.ΩA ≈ m.val.ΩA rtol=1e-12
            @test r.val.Ca ≈ m.val.Ca rtol=1e-12
        end

        for knh3 in ["default", "Clegg"]
            m = whole_system(; measured_kw..., δBT=δBT_MODERN, KNH3_method=knh3)
            r = recalculate_at_target_conditions(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.δBOH₄ ≈ m.val.δBOH₄ rtol=1e-9
            # Looser than the rest, and not because anything is being dropped. A no-op
            # re-solve reuses TA/DIC from the solved state, which have round-tripped
            # through the µmol/mol unit scaling and so differ by ~1 ulp. That reaches H at
            # ~4e-15, and the conversion from fractional abundance to δ has a gain of
            # 1000A/(δ·R_std·(1-A)²) ≈ 271 at A ≈ 0.8, so δBOH₄ lands at ~2e-11. In
            # absolute terms that is 3e-10 ‰, against a measurement precision of ~0.2 ‰.
        end
    end

    @testset "Recalculation matches a direct solve at the target" begin
        # Ground truth: TA and DIC are conservative, so re-solving a measured sample at
        # some other conditions must equal solving there directly with the same settings.
        for (tc, pb) in ((2.0, 0.0), (30.0, 0.0), (20.0, 400.0), (2.0, 400.0))
            for k in ("MyAMI", "Lueker 2000")
                staged = recalculate_at_target_conditions(
                    carbon_system(; measured_kw..., K_method=k); temp_c=tc, pres_bar=pb)
                direct = carbon_system(; measured_kw..., K_method=k, temp_c=tc, pres_bar=pb)
                @test staged.val.pHtot ≈ direct.val.pHtot rtol=1e-12
                @test staged.val.ΩA ≈ direct.val.ΩA rtol=1e-12
            end
        end
    end

    @testset "Settings and inputs are carried on the result" begin
        m = carbon_system(; measured_kw..., K_method="Lueker 2000", KSO4_method="Khoo")
        @test m.settings.K_method == "Lueker 2000"
        @test m.settings.KSO4_method == "Khoo"
        @test m.inputs.temp_c == 20.0
        # A result records the *scope* it was solved with, so a re-solve needs no guesswork
        # about which entry point produced it.
        @test m.system === (:carbon,)
        @test whole_system(; measured_kw..., δBT=δBT_MODERN).system ===
              (:carbon, :boron, :isotopes)

        # ...and survive a recalculation, or the second hop would lose them.
        r = recalculate_at_target_conditions(m; temp_c=2.0)
        @test r.settings.KSO4_method == "Khoo"
        @test r.system === (:carbon,)
    end

    @testset "Inputs the old recursion dropped are now carried" begin
        # Each of these was silently lost when the target state was computed by re-entering
        # the public function with a hand-written argument list.

        # δBT is bulk-conservative: it cannot change with temperature. Uses a non-modern
        # value deliberately - at 39.61 the dropped value and the default coincide, which is
        # exactly how this went unnoticed.
        m = whole_system(; measured_kw..., δBT=38.0)
        r = recalculate_at_target_conditions(m; temp_c=2.0)
        @test r.val.δBT ≈ 38.0
        @test r.val.δBT ≈ m.val.δBT

        # KSO4_method reaches the target conditions.
        khoo = recalculate_at_target_conditions(
            carbon_system(; measured_kw..., K_method="Lueker 2000", KSO4_method="Khoo");
            temp_c=2.0)
        dflt = recalculate_at_target_conditions(
            carbon_system(; measured_kw..., K_method="Lueker 2000");
            temp_c=2.0)
        @test khoo.val.pHsws != dflt.val.pHsws
        @test khoo.val.pHsws ≈ carbon_system(; measured_kw..., temp_c=2.0,
            K_method="Lueker 2000", KSO4_method="Khoo").val.pHsws rtol=1e-12

        # H2ST and NH4T contribute to alkalinity at the target conditions.
        with = recalculate_at_target_conditions(
            carbon_system(; measured_kw..., H2ST=5.0, NH4T=3.0); temp_c=2.0)
        without = recalculate_at_target_conditions(
            carbon_system(; measured_kw...); temp_c=2.0)
        @test with.val.pHtot != without.val.pHtot
    end

    @testset "Only temperature and pressure move" begin
        m = carbon_system(; measured_kw...)
        r = recalculate_at_target_conditions(m; temp_c=2.0, pres_bar=400.0)

        @test r.val.temp_c == 2.0
        @test r.val.pres_bar == 400.0
        @test r.val.sal == m.val.sal          # salinity is not a degree of freedom here
        @test r.val.TA ≈ m.val.TA             # conservative
        @test r.val.DIC ≈ m.val.DIC           # conservative
        @test r.val.pHtot != m.val.pHtot      # ...but the speciation does move

        # An omitted argument keeps its measured value.
        @test recalculate_at_target_conditions(m; temp_c=2.0).val.pres_bar == m.val.pres_bar
        @test recalculate_at_target_conditions(m; pres_bar=400.0).val.temp_c == m.val.temp_c
    end

    @testset "Uncertainties" begin
        # A measurement with uncertainties is a solver built with `varying_errors`; the
        # plain keyword API deliberately has no `errors` argument.
        measured(errs...; kw...) =
            CarbonateSystem(:carbon; varying_errors=errs, measured_kw..., kw...)

        m = measured(:TA, :DIC)(2.0, 2.0)

        # Carried on the result, so not restated.
        r = recalculate_at_target_conditions(m; temp_c=2.0)
        @test r.err.pHtot > 0.0

        # Differentiating the whole chain must agree with solving directly at the target,
        # because TA and DIC are conserved.
        direct = CarbonateSystem(:carbon; varying_errors=(:TA, :DIC),
                                 merge(measured_kw, (temp_c=2.0,))...)(2.0, 2.0)
        @test r.err.pHtot ≈ direct.err.pHtot rtol=1e-8

        # Uncertainty in the collection conditions widens the result.
        @test recalculate_at_target_conditions(m; temp_c=2.0,
                errors=(temp_c=0.5,)).err.pHtot > r.err.pHtot

        # Salinity uncertainty reaches the target state through the constants.
        ms = measured(:TA, :DIC, :sal)(2.0, 2.0, 0.05)
        @test recalculate_at_target_conditions(ms; temp_c=2.0).err.pHtot > r.err.pHtot

        # No uncertainties in, none out.
        @test isnothing(recalculate_at_target_conditions(
            carbon_system(; measured_kw...); temp_c=2.0).err)
    end

    @testset "Presets agree with the solver they wrap" begin
        # `carbon_calculator` / `carbon_boron_calculator` are gone; the four entry points are
        # now presets over CarbonateSystem, and must match calling it directly.
        @test carbon_system(; measured_kw...).val.pHtot ==
              CarbonateSystem(:carbon)(; measured_kw...).val.pHtot
        @test whole_system(; measured_kw..., δBT=δBT_MODERN).val.pHtot ==
              CarbonateSystem(:carbon, :boron, :isotopes)(; measured_kw...,
                                                          δBT=δBT_MODERN).val.pHtot
    end

end
