# Tests for `adjust`: carbon_system/whole_system solve at the conditions they are given, and
# adjust re-solves the same water at another temperature and pressure.

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
                     "Waters 2014", "SB 2020", "Sulpis 2020", "KGen"]

        for k in k_methods
            m = carbon_system(; measured_kw..., K_method=k)
            r = adjust(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.CO₃ ≈ m.val.CO₃ rtol=1e-12
            @test r.val.ΩC ≈ m.val.ΩC rtol=1e-12
        end

        for kso4 in ["default", "Khoo", "WM13"], kf in ["default", "Perez"]
            m = carbon_system(; measured_kw..., K_method="Lueker 2000",
                              KSO4_method=kso4, KF_method=kf)
            r = adjust(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.pHfree ≈ m.val.pHfree rtol=1e-12
        end

        for bt in ["default", "Lee", "KSK18"], ca in ["default", "Culkin", "RT67"]
            m = carbon_system(; measured_kw..., BT_method=bt, Ca_method=ca)
            r = adjust(m)
            @test r.val.pHtot ≈ m.val.pHtot rtol=1e-12
            @test r.val.ΩA ≈ m.val.ΩA rtol=1e-12
            @test r.val.Ca ≈ m.val.Ca rtol=1e-12
        end

        for knh3 in ["default", "Clegg"]
            m = whole_system(; measured_kw..., δBT=δBT_MODERN, KNH3_method=knh3)
            r = adjust(m)
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
            for k in ("KGen", "Lueker 2000")
                staged = adjust(
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
        r = adjust(m; temp_c=2.0)
        @test r.settings.KSO4_method == "Khoo"
        @test r.system === (:carbon,)
    end

    @testset "every input survives a re-solve at target conditions" begin
        # A re-solve inherits its inputs from the original result. Each parameter below is one
        # that a hand-written argument list drops silently, leaving the target state computed
        # from a default instead of the value supplied.

        # δBT is bulk-conservative: it cannot change with temperature. Deliberately non-modern,
        # because at the default 39.61 a dropped value and the default coincide and the test
        # would pass either way.
        m = whole_system(; measured_kw..., δBT=38.0)
        r = adjust(m; temp_c=2.0)
        @test r.val.δBT ≈ 38.0
        @test r.val.δBT ≈ m.val.δBT

        # KSO4_method reaches the target conditions.
        khoo = adjust(
            carbon_system(; measured_kw..., K_method="Lueker 2000", KSO4_method="Khoo");
            temp_c=2.0)
        dflt = adjust(
            carbon_system(; measured_kw..., K_method="Lueker 2000");
            temp_c=2.0)
        @test khoo.val.pHsws != dflt.val.pHsws
        @test khoo.val.pHsws ≈ carbon_system(; measured_kw..., temp_c=2.0,
            K_method="Lueker 2000", KSO4_method="Khoo").val.pHsws rtol=1e-12

        # H2ST and NH4T contribute to alkalinity at the target conditions.
        with = adjust(
            carbon_system(; measured_kw..., H2ST=5.0, NH4T=3.0); temp_c=2.0)
        without = adjust(
            carbon_system(; measured_kw...); temp_c=2.0)
        @test with.val.pHtot != without.val.pHtot
    end

    @testset "Only temperature and pressure move" begin
        m = carbon_system(; measured_kw...)
        r = adjust(m; temp_c=2.0, pres_bar=400.0)

        @test r.val.temp_c == 2.0
        @test r.val.pres_bar == 400.0
        @test r.val.sal == m.val.sal          # salinity is not a degree of freedom here
        @test r.val.TA ≈ m.val.TA             # conservative
        @test r.val.DIC ≈ m.val.DIC           # conservative
        @test r.val.pHtot != m.val.pHtot      # ...but the speciation does move

        # An omitted argument keeps its measured value.
        @test adjust(m; temp_c=2.0).val.pres_bar == m.val.pres_bar
        @test adjust(m; pres_bar=400.0).val.temp_c == m.val.temp_c
    end

    @testset "Uncertainties" begin
        # A measurement with uncertainties is a solver built with `varying_errors`; the
        # plain keyword API deliberately has no `errors` argument.
        measured(errs...; kw...) =
            CarbonateSystem(:carbon; varying_errors=errs, measured_kw..., kw...)

        m = measured(:TA, :DIC)(2.0, 2.0)

        # Carried on the result, so not restated.
        r = adjust(m; temp_c=2.0)
        @test r.err.pHtot > 0.0

        # Differentiating the whole chain must agree with solving directly at the target,
        # because TA and DIC are conserved.
        direct = CarbonateSystem(:carbon; varying_errors=(:TA, :DIC),
                                 merge(measured_kw, (temp_c=2.0,))...)(2.0, 2.0)
        @test r.err.pHtot ≈ direct.err.pHtot rtol=1e-8

        # Uncertainty in the collection conditions widens the result.
        @test adjust(m; temp_c=2.0, σ_temp_c=0.5).err.pHtot > r.err.pHtot

        # Salinity uncertainty reaches the target state through the constants.
        ms = measured(:TA, :DIC, :sal)(2.0, 2.0, 0.05)
        @test adjust(ms; temp_c=2.0).err.pHtot > r.err.pHtot

        # No uncertainties in, none out.
        @test isnothing(adjust(
            carbon_system(; measured_kw...); temp_c=2.0).err)
    end

    @testset "The positional form broadcasts" begin
        # Julia will not broadcast keyword arguments, so the positional form is the only route
        # to a whole cast in one call. It must agree with the keyword form exactly.
        m = carbon_system(; measured_kw...)
        temperatures = [2.0, 6.0, 10.0]
        pressures = [400.0, 200.0, 0.0]

        @test adjust(m, 2.0).pHtot ==
              adjust(m; temp_c=2.0).pHtot
        @test adjust(m, 2.0, 400.0).pHtot ==
              adjust(m; temp_c=2.0, pres_bar=400.0).pHtot

        broadcast_results = adjust.(m, temperatures, pressures)
        loop_results = [adjust(m; temp_c=t, pres_bar=p)
                        for (t, p) in zip(temperatures, pressures)]
        @test length(broadcast_results) == 3
        # Identical, not approximately so: both paths run the same code.
        @test [r.pHtot for r in broadcast_results] == [r.pHtot for r in loop_results]

        # A skipped slot holds that condition at its measured value.
        @test adjust(m, nothing, 400.0).temp_c == m.temp_c

        # A vector of measurements against a vector of collection conditions, which is the
        # case the positional form exists for.
        solve = CarbonateSystem(:carbon; varying=(:TA,), DIC=2000.0, temp_c=20.0, sal=35.0)
        casts = solve.([2300.0, 2310.0, 2290.0])
        @test length(adjust.(casts, temperatures)) == 3
    end

    @testset "Per-row uncertainties broadcast too" begin
        solve = CarbonateSystem(:carbon; varying=(:TA,), varying_errors=(:sal,),
                                DIC=2000.0, temp_c=20.0, sal=35.0)
        casts = solve.([2300.0, 2310.0, 2290.0], 0.05)
        temperatures = [2.0, 6.0, 10.0]
        σ_temperatures = [0.5, 0.2, 0.1]

        results = adjust.(casts, temperatures, nothing, σ_temperatures)

        @test length(results) == 3
        @test all(r -> r.err.pHtot > 0.0, results)
        # Each row carries its own σ, so the propagated uncertainty must differ per row.
        @test length(unique(r.err.pHtot for r in results)) == 3
        # ...and each must match the scalar call it stands for.
        for (cast, result, t, σ) in zip(casts, results, temperatures, σ_temperatures)
            @test result.err.pHtot ==
                  adjust(cast; temp_c=t, σ_temp_c=σ).err.pHtot
        end

        # σ on a condition that was not changed: "the collection temperature equals the
        # measured one, but is uncertain" is legitimate and needs no special case.
        held = adjust(casts[1], nothing, nothing, 0.5)
        @test held.temp_c == 20.0
        @test held.err.pHtot > adjust(casts[1]).err.pHtot
    end

    @testset "Uncertainties in both stages combine in quadrature" begin
        # `temp_c` and `target_temp_c` are separate entries in the AD input vector, so the two
        # thermometers propagate independently. Whether the *measurement* thermometer matters
        # at all depends on which pair was measured, and AD works that out unaided.
        base = (temp_c=25.0, sal=35.0, pres_bar=0.0)

        for (pair, label) in (((TA=2300.0, DIC=2000.0), "TA/DIC"),
                              ((pHtot=8.04, DIC=2000.0), "pH/DIC"))
            inputs = merge(pair, base)
            plain = carbon_system(; inputs...)
            with_measurement_σ = CarbonateSystem(:carbon; varying_errors=(:temp_c,),
                                                 inputs...)(0.1)

            collection_only = adjust(plain; temp_c=2.0, σ_temp_c=0.5).err.pHtot
            measurement_only = adjust(with_measurement_σ; temp_c=2.0).err.pHtot
            both = adjust(with_measurement_σ;
                          temp_c=2.0, σ_temp_c=0.5).err.pHtot

            @test both ≈ sqrt(collection_only^2 + measurement_only^2) rtol=1e-9
        end

        # With TA and DIC given, they are conservative and carried across unchanged, so the
        # measurement temperature never reaches the collection state. This is the sentinel: if
        # it ever becomes materially nonzero, something is leaking out of stage 1.
        conservative = CarbonateSystem(:carbon; varying_errors=(:temp_c,),
                                       TA=2300.0, DIC=2000.0, base...)(0.1)
        @test adjust(conservative; temp_c=2.0).err.pHtot < 1e-10

        # With pH measured, TA is *derived* at the measurement temperature, so it does.
        derived = CarbonateSystem(:carbon; varying_errors=(:temp_c,),
                                  pHtot=8.04, DIC=2000.0, base...)(0.1)
        @test adjust(derived; temp_c=2.0).err.pHtot > 1e-4
    end

    @testset "Arguments adjust does not change are refused" begin
        m = carbon_system(; measured_kw...)

        # Salinity: the mistake people will actually make, so it explains itself rather than
        # reporting an unknown name.
        @test_throws "different water" adjust(m; sal=34.0)
        @test_throws "direct solve" adjust(m; sal=34.0)
        @test_throws "varying_errors" adjust(m; σ_sal=0.05)

        # `errors=` is what this function used to take, so it names its replacement.
        @test_throws "σ_temp_c and σ_pres_bar" adjust(m;
                                                                       errors=(temp_c=0.5,))

        @test_throws ArgumentError adjust(m; temperature=2.0)

        # A vector handed to the keyword form used to die as `isless(::Int64, ::Vector)` from
        # inside the solver, or `convert(Float64, ::Vector)` from the uncertainty path.
        @test_throws "Broadcast the positional form" adjust(m;
                                                                             temp_c=[2.0, 4.0])
        @test_throws "Broadcast the positional form" adjust(m;
                                                                             σ_temp_c=[0.5, 0.2])
        @test_throws "Broadcast the positional form" adjust(m, [2.0, 4.0])
    end

    @testset "An adjusted result cannot be adjusted again" begin
        m = carbon_system(; measured_kw...)
        collected = adjust(m; temp_c=2.0)

        # Having a `source` is what marks a result as already carried; there is no separate flag.
        @test isnothing(m.source)
        @test !isnothing(collected.source)

        @test_throws "produced by `adjust`" adjust(
            collected; temp_c=4.0)
        # Positional form and no-argument form go through the same guard.
        @test_throws ArgumentError adjust(collected, 4.0)
        @test_throws ArgumentError adjust(collected)

        # Carrying the *measurement* somewhere else is the supported route, and stays open.
        @test adjust(m; temp_c=4.0).pHtot != collected.pHtot

        # The source is recorded down the uncertainty path too, which builds the result on a
        # different branch — and the intermediate measured state inside it must not get one.
        uncertain = CarbonateSystem(:carbon; varying_errors=(:TA,), measured_kw...)(2.0)
        @test !isnothing(adjust(uncertain; temp_c=2.0, σ_temp_c=0.5).source)
        @test_throws ArgumentError adjust(
            adjust(uncertain; temp_c=2.0); temp_c=4.0)
    end

    @testset "Target conditions are validated" begin
        # Both exits of stage 2 call the core directly rather than going through `_run`, so
        # before this check existed a 200 °C target or a negative pressure returned a number
        # in silence while stage 1 warned or threw for the same values.
        m = carbon_system(; measured_kw...)

        @test_throws ArgumentError adjust(m; pres_bar=-50.0)
        @test_throws ArgumentError adjust(m; temp_c=NaN)

        # Out of range warns rather than throwing: the package should still compute for a
        # hydrothermal vent, so the value has to come back as well as the warning.
        hot = @test_logs (:warn,) match_mode=:any adjust(m;
                                                                                  temp_c=200.0)
        @test hot.temp_c == 200.0

        # In range stays quiet, including on the uncertainty path — the check runs before the
        # AD call, where the values are still plain numbers rather than Duals.
        @test_logs adjust(m; temp_c=2.0)
        uncertain = CarbonateSystem(:carbon; varying_errors=(:TA,), measured_kw...)(2.0)
        @test_logs adjust(uncertain; temp_c=2.0)

        # Known blind spot, pinned rather than fixed: `_check_conditions` returns early for
        # anything that is not `Real`, and `missing isa Real` is false, so a `missing` target
        # still fails from deep inside the solver instead of being named here. The exception
        # type depends on which parameterisation the solve reaches, so only the failure is
        # asserted, not its shape.
        @test_throws Exception adjust(m; temp_c=missing)
    end

    @testset "Only a non-conservative pair repeats the original solve" begin
        # From TA and DIC, nothing at the original conditions can reach the adjusted state,
        # so uncertainty in the original temperature adds nothing.
        from_TA_DIC = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC, :temp_c),
                                      measured_kw...)(2.0, 2.0, 0.5)
        without_T = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC),
                                    measured_kw...)(2.0, 2.0)
        @test adjust(from_TA_DIC; temp_c = 2.0).err.pHtot ≈
              adjust(without_T; temp_c = 2.0).err.pHtot rtol = 1e-10

        # From pH and TA, DIC was derived with constants at the original temperature, so its
        # uncertainty does reach the adjusted state.
        ph_kw = (TA = 2300.0, pHtot = 8.0, temp_c = 20.0, sal = 35.0, pres_bar = 0.0)
        from_pH = CarbonateSystem(:carbon; varying_errors = (:TA, :pHtot, :temp_c),
                                  ph_kw...)(2.0, 0.005, 0.5)
        without_T_pH = CarbonateSystem(:carbon; varying_errors = (:TA, :pHtot),
                                       ph_kw...)(2.0, 0.005)
        @test adjust(from_pH; temp_c = 2.0).err.pHtot > adjust(without_T_pH; temp_c = 2.0).err.pHtot

        # Without uncertainties, adjusting matches solving directly at the new conditions.
        @test adjust(carbon_system(; measured_kw...); temp_c = 4.0, pres_bar = 150.0).pHtot ≈
              carbon_system(; measured_kw..., temp_c = 4.0, pres_bar = 150.0).pHtot rtol = 1e-12
    end

    @testset "A result without carbon cannot be adjusted" begin
        @test_throws ArgumentError adjust(boron_system(pHtot = 8.1); temp_c = 2.0)
    end

    @testset "The former name still works" begin
        m = carbon_system(; measured_kw...)
        @test (@test_deprecated at_collection_conditions(m; temp_c = 2.0)).pHtot ==
              adjust(m; temp_c = 2.0).pHtot
        @test (@test_deprecated at_collection_conditions(m, 2.0, 400.0)).pHtot ==
              adjust(m, 2.0, 400.0).pHtot
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
