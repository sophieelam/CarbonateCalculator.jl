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
