# Choosing the constants with `K_settings`, and the warnings that go with it.
#
# Three kinds of warning, each tested both ways — fires when it should, silent when it should
# not, because a warning that fires on correct input trains people to ignore it:
#
#   companions   a K1/K2 parameterisation selected without the constants it was published
#                with, a companion-only option used with another, or a choice Kgen ignores.
#                Depend only on the settings, so they fire when the settings are made.
#   range        the K1/K2 parameterisation used outside the temperature and salinity it was
#                fitted over, or with a seawater composition only KGen corrects for. Depend
#                on the water, so they fire per calculation.
#   legacy       the old flat keywords (`K_method = …`), which still work and say what to
#                write instead.
#
# And the property that motivated the change: no choice silently changes another.

using Test
using Logging
using CarbonateCalculator
import CarbonateCalculator: calculate_constants
import CarbonateCalculator.Constants: Takahashi_fH, Uppstrom_BT, case67_BT, MODERN_CALCIUM,
                                      K1_K2_RANGES, SPECIATION

const SEAWATER = (TA = 2300.0, DIC = 2000.0, temp_c = 20.0, sal = 35.0)
const FRESHWATER_SET = (K1_K2_method = "Millero 1979", fH_method = "freshwater",
                        BT_method = "freshwater", pressure_method = "Millero 1983")
const GEOSECS_SET = (K1_K2_method = "Mehrbach 1973", BT_method = "Culkin",
                     Ca_method = "Culkin")
const PENG_SET = (K1_K2_method = "Mehrbach 1973", fH_method = "Peng 1987",
                  BT_method = "Culkin", Ca_method = "Culkin")

"Build a KSettings with no warnings allowed."
silent_settings(; kw...) = @test_logs min_level=Logging.Warn KSettings(; kw...)

@testset "K_settings" begin

    @testset "building settings" begin
        default = KSettings()
        @test default.K1_K2_method == "KGen"
        @test default.KSO4_method == "Dickson"
        @test default.KF_method == "Dickson"
        @test default.KNH3_method == "Millero"
        @test default.BT_method == "Uppstrom"
        @test default.Ca_method == "modern"
        @test default.fH_method == "Takahashi 1982"
        @test default.pressure_method == "Millero 1995"
        @test default.MyAMI_mode == "approximate"

        # Only what differs from the defaults is named; the rest keep theirs.
        roy = KSettings(K1_K2_method = "Roy 1993", KF_method = "Perez")
        @test roy.K1_K2_method == "Roy 1993" && roy.KF_method == "Perez"
        @test roy.KSO4_method == "Dickson"

        # "default" and Symbols are accepted, and resolve to the same thing.
        @test KSettings(K1_K2_method = :default, KF_method = :Perez) ==
              KSettings(KF_method = "Perez")
        @test KSettings(KSO4_method = "default") == default

        # Copy with changes.
        @test KSettings(roy; KF_method = "Dickson") == KSettings(K1_K2_method = "Roy 1993")
        @test KSettings(roy; KNH3_method = "Clegg").KF_method == "Perez"

        # Printed as the call that rebuilds it, naming only the changes.
        @test repr(default) == "KSettings()"
        @test repr(roy) == "KSettings(K1_K2_method = \"Roy 1993\", KF_method = \"Perez\")"
        @test eval(Meta.parse(repr(roy))) == roy

        # Equal settings hash equally, so they work as Dict keys and warning ids.
        @test hash(KSettings(KF_method = "Perez")) == hash(KSettings(KF_method = :Perez))
    end

    @testset "bad settings are errors, naming the fix" begin
        function err(f)
            try
                f()
            catch e
                e isa ArgumentError && return e.msg
                rethrow()
            end
            return ""
        end

        @test_throws ArgumentError KSettings(K1_K2_method = "Lueker 2001")
        @test occursin("Lueker 2000", err(() -> KSettings(K1_K2_method = "Lueker 2001")))
        @test_throws ArgumentError KSettings(KF_method = "Dixon")
        @test_throws ArgumentError KSettings(BT_method = 3)

        # Old K_method is renamed, and a near-miss field name gets a suggestion.
        @test occursin("K1_K2_method", err(() -> KSettings(K_method = "Roy 1993")))
        @test occursin("did you mean KF_method", err(() -> KSettings(kf_method = "Perez")))
        @test_throws ArgumentError KSettings(K1K2 = "Roy 1993")

        # The merged and removed names say what replaced them.
        msg = err(() -> KSettings(K1_K2_method = "Mehrbach 1973 A"))
        @test occursin("\"Mehrbach 1973\"", msg) && occursin("Culkin", msg)
        msg = err(() -> KSettings(K1_K2_method = "Mehrbach 1973 B"))
        @test occursin("Peng 1987", msg)
        @test occursin("README", err(() -> KSettings(K1_K2_method = "automatic")))

        # Only the approximation to MyAMI exists in Kgen.jl.
        @test occursin("approximation", err(() -> KSettings(MyAMI_mode = "calculate")))

        # Wrong type for K_settings itself.
        @test_throws ArgumentError carbon_system(; SEAWATER..., K_settings = "Roy 1993")
        @test_throws ArgumentError CarbonateSystem(:carbon; varying = (:TA, :DIC),
                                                   K_settings = (K1_K2_method = "x",))
    end

    @testset "companion warnings" begin
        # Defaults and every self-contained parameterisation are silent.
        silent_settings()
        for (name, row) in SPECIATION.K1_K2_method.options
            isempty(row.companions) && silent_settings(K1_K2_method = name)
        end

        # Mehrbach 1973 alone names both published sets and what is selected instead.
        @test_logs (:warn, r"You've chosen \"Mehrbach 1973\" for K1_K2_method.*intended to be used with one of these sets.*GEOSECS.*Peng.*You currently have \"Takahashi 1982\" selected for fH_method, \"Uppstrom\" selected for BT_method, \"modern\" selected for Ca_method"s) KSettings(K1_K2_method = "Mehrbach 1973")
        # A partial set is still not the set.
        @test_logs (:warn, r"Mehrbach 1973.*\"Culkin\" selected for BT_method, \"modern\" selected for Ca_method"s) KSettings(K1_K2_method = "Mehrbach 1973", BT_method = "Culkin")
        # Either complete set is silent.
        silent_settings(; GEOSECS_SET...)
        silent_settings(; PENG_SET...)

        # Millero 1979 likewise, with its single freshwater set.
        @test_logs (:warn, r"\"Millero 1979\".*intended to be used with this set.*freshwater"s) KSettings(K1_K2_method = "Millero 1979")
        silent_settings(; FRESHWATER_SET...)
        # Extra choices beyond the set do not matter to it.
        silent_settings(; FRESHWATER_SET..., KNH3_method = "Clegg")

        # Companion-only options used with any other K1/K2 parameterisation.
        @test_logs (:warn, r"chosen \"Peng 1987\" for fH_method.*intended to be used with K1_K2_method = \"Mehrbach 1973\".*currently have \"Lueker 2000\""s) KSettings(K1_K2_method = "Lueker 2000", fH_method = "Peng 1987")
        @test_logs (:warn, r"chosen \"freshwater\" for BT_method.*\"Millero 1979\""s) KSettings(BT_method = "freshwater", K1_K2_method = "Roy 1993")
        @test_logs (:warn, r"chosen \"Millero 1983\" for pressure_method"s) KSettings(pressure_method = "Millero 1983", K1_K2_method = "Roy 1993")
        # Culkin boron is a parameterisation in its own right, and needs no particular K1/K2.
        silent_settings(BT_method = "Culkin", K1_K2_method = "Roy 1993")

        # KGen computes KSO4, KF and its own pressure correction, so choosing them is a no-op.
        @test_logs (:warn, r"chosen \"Perez\" for KF_method, but .* \"KGen\" selected .* Kgen computes KF itself"s) KSettings(KF_method = "Perez")
        @test_logs (:warn, r"chosen \"Khoo\" for KSO4_method.*Kgen computes KSO4 itself"s) KSettings(KSO4_method = "Khoo")
        @test_logs (:warn, r"only reaches KH2S and KNH3"s) match_mode=:any KSettings(pressure_method = "Millero 1983")
        # ...but not for any other K1/K2 parameterisation, where they are used.
        silent_settings(K1_K2_method = "Lueker 2000", KF_method = "Perez", KSO4_method = "Khoo")
        # KNH3 and BT are used under KGen too.
        silent_settings(KNH3_method = "Clegg", BT_method = "Lee")

        # The same warnings reach the solver, once, when it is built — not per sample.
        solve = @test_logs (:warn, r"Mehrbach 1973") CarbonateSystem(:carbon;
            varying = (:TA, :DIC), K_settings = (K1_K2_method = "Mehrbach 1973",))
        @test_logs min_level=Logging.Warn solve.([2300.0, 2310.0], [2000.0, 2010.0])
        # ...and through a keyword call, which builds its settings per call.
        @test_logs (:warn, r"Mehrbach 1973") match_mode=:any carbon_system(; SEAWATER...,
            K_settings = (K1_K2_method = "Mehrbach 1973",))
    end

    @testset "no choice changes another" begin
        # The point of the change. Mehrbach 1973 used to bring Culkin boron with it, and
        # Millero 1979 a whole freshwater environment; now each takes the defaults unless told.
        default = carbon_system(; SEAWATER...)
        mehrbach = carbon_system(; SEAWATER..., K_settings = (K1_K2_method = "Mehrbach 1973",))
        @test mehrbach.BT == default.BT
        @test mehrbach.BT ≈ 1e6 * Uppstrom_BT(sal = 35.0).BT
        geosecs = carbon_system(; SEAWATER..., K_settings = GEOSECS_SET)
        @test geosecs.BT ≈ 1e6 * case67_BT(sal = 35.0).BT
        @test geosecs.BT != default.BT

        bare = calculate_constants(temp_c = 20.0, sal = 0.0,
                                   K_settings = (K1_K2_method = "Millero 1979",))
        @test bare.fH == Takahashi_fH(temp_c = 20.0, sal = 0.0).fH
        fresh = calculate_constants(temp_c = 20.0, sal = 0.0, K_settings = FRESHWATER_SET)
        @test fresh.fH == 1.0
        @test fresh.BT == 0.0
        # The K1/K2 do not depend on the companions at the surface, only on fH...
        @test fresh.Ks.K1 != bare.Ks.K1
        # ...and the pressure table is its own choice too.
        deep(pressure_method) = calculate_constants(temp_c = 20.0, sal = 0.0, pres_bar = 300.0,
            K_settings = merge(FRESHWATER_SET, (; pressure_method))).Ks.K1
        @test deep("Millero 1983") != deep("Millero 1995")

        # The settings used are reported, as given.
        @test fresh.method == KSettings(; FRESHWATER_SET...)
        @test mehrbach.settings.K_settings == KSettings(K1_K2_method = "Mehrbach 1973")
    end

    @testset "range warnings" begin
        lueker = (K1_K2_method = "Lueker 2000",)

        # Inside the fitted range: silent, at both edges.
        for temp_c in (2.0, 20.0, 35.0), sal in (19.0, 35.0, 43.0)
            @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., temp_c, sal,
                                                             K_settings = lueker)
        end

        # Outside it: the user's wording, the value, the range and where it comes from.
        @test_logs (:warn, r"You currently have selected \"Lueker 2000\" for K1_K2_method, however this method is not recommended for temp_c = 40\.0 °C: it was fitted over 2\.0 to 35\.0 °C \(Lueker et al\. \(2000\)"s) carbon_system(; SEAWATER..., temp_c = 40.0, K_settings = lueker)
        @test_logs (:warn, r"not recommended for sal = 10\.0: it was fitted over 19\.0 to 43\.0"s) carbon_system(; SEAWATER..., sal = 10.0, K_settings = lueker)
        @test_logs (:warn, r"not recommended for temp_c = 0\.5"s) carbon_system(; SEAWATER..., temp_c = 0.5, K_settings = lueker)

        # Every parameterisation has a range, and every range is checked.
        @test Set(keys(K1_K2_RANGES)) == Set(keys(SPECIATION.K1_K2_method.options))
        for (name, r) in K1_K2_RANGES
            @test r.temp_c[1] <= r.temp_c[2] && r.sal[1] <= r.sal[2]
            @test !isempty(r.source)
        end

        # Millero 1979 is pure water only.
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., sal = 0.0,
                                                         K_settings = FRESHWATER_SET)
        @test_logs (:warn, r"\"Millero 1979\".*not recommended for sal = 35\.0: it was fitted over sal = 0\.0 only"s) carbon_system(; SEAWATER..., K_settings = FRESHWATER_SET)

        # A brine Papadimitriou 2018 covers raises nothing at all — neither its own range
        # check nor the global one, which is the union of the registry's ranges...
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., temp_c = 10.0,
            sal = 60.0, K_settings = (K1_K2_method = "Papadimitriou 2018",))
        # ...while the default, which is not fitted there, says so.
        @test_logs (:warn, r"\"KGen\".*sal = 60\.0"s) carbon_system(; SEAWATER...,
                                                                  temp_c = 10.0, sal = 60.0)

        # The default, KGen, has a range too.
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER...)
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., temp_c = 0.0)
        @test_logs (:warn, r"\"KGen\".*temp_c = -3\.0"s) carbon_system(; SEAWATER...,
                                                                     temp_c = -3.0)

        # Once per method and side, not once per sample.
        solve = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c), K_settings = lueker)
        @test_logs (:warn, r"temp_c = 40\.0") solve.([2300.0, 2300.0, 2300.0],
                                                    [2000.0, 2000.0, 2000.0],
                                                    [40.0, 41.0, 42.0])
        @test_logs (:warn, r"temp_c = 40\.0") (:warn, r"temp_c = 0\.0") solve.(
            [2300.0, 2300.0, 2300.0], [2000.0, 2000.0, 2000.0], [40.0, 0.0, 41.0])

        # adjust checks the new conditions, which is where an extrapolation usually happens.
        measured = carbon_system(; SEAWATER..., K_settings = lueker)
        @test_logs (:warn, r"temp_c = 0\.5") adjust(measured; temp_c = 0.5)
        @test_logs min_level=Logging.Warn adjust(measured; temp_c = 5.0)

        # Constants supplied ready-made were not computed from the settings, so the settings
        # are not checked against them.
        env = calculate_constants(temp_c = 40.0, sal = 35.0, K_settings = lueker)
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., temp_c = 40.0,
                                                         Ks = env, K_settings = lueker)
    end

    @testset "composition warnings" begin
        # Only KGen corrects the constants for Ca and Mg.
        @test_logs (:warn, r"\"Lueker 2000\".*not recommended for a non-modern seawater composition.*only \"KGen\" corrects"s) carbon_system(; SEAWATER..., Ca = 0.02, K_settings = (K1_K2_method = "Lueker 2000",))
        @test_logs (:warn, r"non-modern"s) carbon_system(; SEAWATER..., Mg = 0.03, K_settings = (K1_K2_method = "Roy 1993",))
        # Close to modern, so this is about the composition check and not Kgen's own range.
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., Ca = 0.011)
        # An explicitly modern composition is not non-modern.
        @test_logs min_level=Logging.Warn carbon_system(; SEAWATER..., Ca = MODERN_CALCIUM,
            K_settings = (K1_K2_method = "Lueker 2000",))
    end

    @testset "the old keywords still work, and say what to write" begin
        lueker = carbon_system(; SEAWATER..., K_settings = (K1_K2_method = "Lueker 2000",))
        old = @test_logs (:warn, r"K_method is replaced by K_settings.*K_settings = \(K1_K2_method = \"Lueker 2000\",\)"s) carbon_system(; SEAWATER..., K_method = "Lueker 2000")
        @test old.pHtot == lueker.pHtot
        @test old.settings.K_settings == lueker.settings.K_settings

        # The three that used to switch companions silently translate to them explicitly, so
        # an old script gets the same numbers — and is told which settings those were.
        translated = (
            "Mehrbach 1973 A" => (K1_K2_method = "Mehrbach 1973", BT_method = "Culkin"),
            "Mehrbach 1973 B" => (K1_K2_method = "Mehrbach 1973", fH_method = "Peng 1987"),
            "Millero 1979" => FRESHWATER_SET,
            "default" => (K1_K2_method = "KGen",),
        )
        for (name, settings) in translated
            water = name == "Millero 1979" ? merge(SEAWATER, (sal = 0.0,)) : SEAWATER
            old = with_logger(NullLogger()) do
                carbon_system(; water..., K_method = name)
            end
            new = with_logger(NullLogger()) do
                carbon_system(; water..., K_settings = settings)
            end
            @test old.pHtot == new.pHtot
            @test old.BT == new.BT
            @test old.settings.K_settings == KSettings(; settings...)
        end

        # The translation of K_method wins over an old BT_method, as the override used to.
        mixed = with_logger(NullLogger()) do
            carbon_system(; SEAWATER..., K_method = "Mehrbach 1973 A", BT_method = "Lee")
        end
        @test mixed.settings.K_settings.BT_method == "Culkin"

        # Several old keywords at once, including on a solver.
        solve = @test_logs (:warn, r"K_method, KF_method are replaced"s) match_mode=:any CarbonateSystem(:carbon; varying = (:TA, :DIC), K_method = "Roy 1993", KF_method = "Perez")
        @test solve.settings.K_settings == KSettings(K1_K2_method = "Roy 1993", KF_method = "Perez")

        # calculate_constants takes them too.
        env = @test_logs (:warn, r"replaced by K_settings") calculate_constants(
            temp_c = 20.0, sal = 35.0, K_method = "Lueker 2000")
        @test env.Ks == calculate_constants(temp_c = 20.0, sal = 35.0,
                                            K_settings = (K1_K2_method = "Lueker 2000",)).Ks
        @test_throws ArgumentError calculate_constants(temp_c = 20.0, sal = 35.0, K_methodd = "x")

        # Errors: "automatic" is gone, both forms at once is ambiguous, and settings cannot
        # vary within one solver.
        @test_throws ArgumentError carbon_system(; SEAWATER..., K_method = "automatic")
        @test_throws ArgumentError carbon_system(; SEAWATER..., K_method = "Roy 1993",
                                                 K_settings = (KF_method = "Perez",))
        @test_throws ArgumentError CarbonateSystem(:carbon; varying = (:TA, :DIC, :K_method))
    end

    @testset "one parameterisation across the domain" begin
        # The default does not switch parameterisation with the conditions, so pH is smooth
        # in temperature — the step `"automatic"` used to put at 2 °C is gone with it.
        pH(t) = with_logger(NullLogger()) do
            carbon_system(; SEAWATER..., temp_c = t).pHtot
        end
        @test abs(pH(2.0) - pH(1.99)) < 1e-3
    end
end
