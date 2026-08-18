# Broadcasting a solver.
#
# Julia cannot broadcast over keyword arguments, so processing many samples used to mean a
# comprehension. Naming the varying parameters up front hands back a scalar function, which
# broadcasts — and checks the names when it is built rather than on every row.
#
# There is no `carbon_solver`: a solver is `CarbonateSystem(scope; varying = …)`, which reads
# the same way for every scope.

using Test
using CarbonateCalculator

@testset "Solvers" begin

    TA   = [2300.0, 2310.0, 2290.0]
    DIC  = [2000.0, 2010.0, 1990.0]
    temp = [25.0, 10.0, 2.0]
    sal  = [35.0, 34.5, 34.0]

    @testset "broadcasting agrees with the equivalent loop" begin
        solve = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c, :sal),
                                K_method = "Lueker 2000")
        broadcast_results = solve.(TA, DIC, temp, sal)
        loop_results = [carbon_system(TA = t, DIC = d, temp_c = tc, sal = s,
                                      K_method = "Lueker 2000")
                        for (t, d, tc, s) in zip(TA, DIC, temp, sal)]

        @test length(broadcast_results) == 3
        for (a, b) in zip(broadcast_results, loop_results)
            # Identical, not approximately so: both paths run the same code.
            @test a.pHtot === b.pHtot
            @test a.CO₃ === b.CO₃
            @test a.revelle_factor === b.revelle_factor
        end
    end

    @testset "results stay concretely typed" begin
        # If the element type were abstract, every downstream access would be a dynamic
        # dispatch and vectorised work would be pointless.
        results = CarbonateSystem(:carbon; varying = (:TA, :DIC)).(TA, DIC)
        @test isconcretetype(eltype(results))
        @test Base.return_types(rs -> [r.pHtot for r in rs], (typeof(results),))[1] ===
              Vector{Float64}
    end

    @testset "scalars and settings" begin
        solve = CarbonateSystem(:carbon; varying = (:TA, :DIC), temp_c = 2.0, sal = 34.5)
        # A fixed setting applies to every element and cannot drift between rows.
        @test solve(2300.0, 2000.0).temp_c == 2.0
        @test all(r.sal == 34.5 for r in solve.(TA, DIC))
        # Scalars broadcast against vectors as usual.
        @test length(solve.(TA, 2000.0)) == 3
    end

    @testset "names are checked when the solver is built" begin
        # The value of doing this up front is that a mistake surfaces immediately rather
        # than on row 700,000 of a broadcast.
        C(; kw...) = CarbonateSystem(:carbon; kw...)
        @test_throws ArgumentError C(varying = (:TA, :DICC))           # misspelled
        @test_throws ArgumentError C(varying = (:TA,))                 # under-determined
        @test_throws ArgumentError C(varying = (:TA, :DIC, :pHtot))    # over-determined
        @test_throws ArgumentError C(varying = (:TA, :pHtot, :pHsws))  # one freedom twice
        @test_throws ArgumentError C(varying = (:TA, :DIC, :TA))       # named twice
        @test_throws ArgumentError C(varying = (:TA, :DIC), TA = 2300.0)  # varying and fixed
        @test_throws ArgumentError C(varying = (:TA, :DIC), K_methodd = "x")  # bad setting
        # ...and a parameter the scope does not cover.
        @test_throws ArgumentError C(varying = (:TA, :DIC), δBT = 39.61)
    end

    @testset "calling it wrongly" begin
        solve = CarbonateSystem(:carbon; varying = (:TA, :DIC))
        @test_throws ArgumentError solve(2300.0)
        @test_throws ArgumentError solve(2300.0, 2000.0, 25.0)
        # A solver with no varying parameters is the keyword form, and says so.
        @test_throws ArgumentError CarbonateSystem(:carbon)(2300.0, 2000.0)
    end

    @testset "uncertainties" begin
        # Uncertainties are named at construction and supplied positionally, after the
        # varying parameters. Every one of them varies: a scalar broadcasts against a
        # vector, so a σ that is the same for every sample needs no separate mechanism.
        solve = CarbonateSystem(:carbon; varying = (:TA, :DIC),
                                varying_errors = (:TA, :DIC))
        results = solve.(TA, DIC, 2.0, 2.0)
        @test all(!isnothing(r.err) for r in results)
        @test all(r.err.pHtot > 0 for r in results)
        # ...and the uncertainties must not leak into the recorded inputs.
        @test !haskey(results[1].inputs, :errors)
        @test !haskey(results[1].inputs, :TA_error)

        # A per-sample σ column gives a different answer per row, where a scalar does not.
        per_sample = solve.(TA, DIC, 2.0, [1.0, 5.0, 20.0])
        @test per_sample[1].err.pHtot < per_sample[3].err.pHtot

        # The plain keyword API has no `errors` at all — propagation is reached by building
        # a solver, not by passing a NamedTuple to `carbon_system`.
        @test_throws ArgumentError carbon_system(TA = 2300.0, DIC = 2000.0,
                                                 errors = (TA = 2.0,))

        # A misspelled uncertainty fails when the solver is built.
        @test_throws ArgumentError CarbonateSystem(:carbon; varying = (:TA, :DIC),
                                                   varying_errors = (:TAA,))
        # ...as does one naming a parameter this solver has no value for.
        @test_throws ArgumentError CarbonateSystem(:carbon; varying = (:TA, :DIC),
                                                   varying_errors = (:pCO₂,))
        # A σ for a fixed setting is fine: the value exists, it just does not vary.
        @test CarbonateSystem(:carbon; varying = (:TA, :DIC), temp_c = 20.0,
                              varying_errors = (:temp_c,))(2300.0, 2000.0, 0.5).err.pHtot > 0
    end

    @testset "every input pair propagates a condition uncertainty" begin
        # An uncertainty on temp_c reaches the answer only through the equilibrium
        # constants, so it exercises a part of each solve that a σ on TA or DIC does not.
        # Both ways that can break are silent rather than loud:
        #
        #   - a bracketing `find_zero` returns a Float64 whatever it is given, so the
        #     ForwardDiff partials are dropped and σ comes back as exactly 0.0;
        #   - a Newton solve whose initial guess is typed from the concentrations alone
        #     meets a Dual residual with a Float64 state and throws from inside Roots.
        #
        # Neither was covered, because every uncertainty test carried a σ on TA or DIC too,
        # which made the solver state Dual for the wrong reason. Asserting on temp_c alone
        # is what distinguishes the two.
        pairs = [(:TA, :DIC)    => (2300.0, 2000.0),
                 (:CO₂, :DIC)   => (10.0, 2000.0),
                 (:CO₃, :DIC)   => (200.0, 2000.0),
                 (:HCO₃, :CO₃)  => (1800.0, 200.0),
                 (:HCO₃, :DIC)  => (1800.0, 2000.0),
                 (:CO₂, :HCO₃)  => (10.0, 1800.0),
                 (:CO₂, :CO₃)   => (10.0, 200.0),
                 (:CO₂, :TA)    => (10.0, 2300.0),
                 (:HCO₃, :TA)   => (1800.0, 2300.0),
                 (:CO₃, :TA)    => (200.0, 2300.0)]

        for (names, values) in pairs
            solve = CarbonateSystem(:carbon; temp_c = 20.0, varying_errors = (:temp_c,),
                                    NamedTuple{names}(values)...)
            σ = solve(0.5).err.pHtot
            @test σ > 0            # 0.0 means the partials were silently discarded
            @test isfinite(σ)
        end
    end

    @testset "partial uncertainties combine in quadrature" begin
        # A parameter with no entry in `varying_errors` is treated as exact, contributing
        # zero variance — the standard first-order result, and worth pinning because it is
        # easy to assume an unlisted parameter is being accounted for.
        with(errs...) = CarbonateSystem(:carbon; varying = (:TA, :DIC), temp_c = 20.0,
                                        varying_errors = errs)

        from_DIC  = with(:DIC)(2300.0, 2000.0, 2.0).err.pHtot
        from_TA   = with(:TA)(2300.0, 2000.0, 2.0).err.pHtot
        from_both = with(:DIC, :TA)(2300.0, 2000.0, 2.0, 2.0).err.pHtot
        @test sqrt(from_DIC^2 + from_TA^2) ≈ from_both rtol=1e-9

        # Adding an uncertainty can only widen the result.
        with_temp = with(:DIC, :TA, :temp_c)(2300.0, 2000.0, 2.0, 2.0, 0.5).err.pHtot
        @test with_temp > from_both
    end

    @testset "every scope can be a solver" begin
        # The reason there is no carbon_solver/whole_solver pair: this reads the same for all
        # four scopes, where those covered only two.
        whole = CarbonateSystem(:carbon, :boron, :isotopes;
                                varying = (:TA, :DIC), δBT = 39.61)
        @test whole(2300.0, 2000.0).pHtot ≈
              whole_system(TA = 2300.0, DIC = 2000.0, δBT = 39.61).pHtot

        boron = CarbonateSystem(:boron, :isotopes; varying = (:BT, :BOH₄), δBT = 39.61)
        @test boron(415.7, 100.0).pHtot ≈
              boron_system(BT = 415.7, BOH₄ = 100.0, δBT = 39.61).pHtot

        isotopes = CarbonateSystem(:isotopes; varying = (:pHtot,), δBT = 39.61)
        @test isotopes(8.1).δBOH₄ ≈ boron_isotopes(pHtot = 8.1, δBT = 39.61).δBOH₄
    end

end
