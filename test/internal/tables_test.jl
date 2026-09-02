# A vector of results is a Tables.jl table, one row per result.
#
# Before this interface existed, `DataFrame(results)` did not fail — it fell back to reflecting
# over `propertynames`, which advertises the struct fields as well as the computed values, and
# produced a 47-column frame carrying `val` twice over (once nested, once flattened) beside
# columns of `Ks`, `inputs` and `nothing`. The assertions below are what replaced that.

using Test
using CarbonateCalculator
using DataFrames
using Tables

@testset "Tables.jl interop" begin

    temperatures = [10.0, 20.0]
    plain = [carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = t) for t in temperatures]

    solve = CarbonateSystem(:carbon; varying = (:temp_c,), varying_errors = (:TA, :DIC),
                            TA = 2300.0, DIC = 2000.0)
    uncertain = solve.(temperatures, 2.0, 2.0)

    @testset "a vector of results declares itself a table" begin
        @test Tables.istable(typeof(plain))
        @test Tables.rowaccess(typeof(plain))
        # A broadcast gives a concrete element type, which is what lets the ragged-row filling
        # below be skipped entirely on the path that matters.
        @test isconcretetype(eltype(uncertain))
        # The schema belongs to the rows, not to the vector.
        @test Tables.schema(Tables.rows(plain)).names == keys(first(plain).val)
    end

    @testset "without uncertainties the row is val" begin
        df = DataFrame(plain)

        @test ncol(df) == length(keys(first(plain).val))
        @test df == DataFrame([r.val for r in plain])   # the idiom this replaces
        @test all(t -> t === Float64, eltype.(eachcol(df)))
        @test df.pHtot == [r.pHtot for r in plain]

        # The struct fields are gone: they were never columns, only tab-completion entries.
        @test !("val" in names(df))
        @test !("err" in names(df))
        @test !("Ks" in names(df))
    end

    @testset "with uncertainties each value gains a σ_ column" begin
        df = DataFrame(uncertain)
        computed = keys(first(uncertain).val)

        @test ncol(df) == 2 * length(computed)
        @test all(name -> "σ_$name" in names(df), computed)
        @test all(t -> t === Float64, eltype.(eachcol(df)))

        @test df.pHtot == [r.pHtot for r in uncertain]
        @test df.σ_pHtot == [r.err.pHtot for r in uncertain]
        # Uncertainties vary per row here, so a constant column would pass the check above
        # without the σ actually tracking anything.
        @test df.σ_pHtot[1] != df.σ_pHtot[2]
    end

    @testset "a row missing an uncertainty gets nothing, not a dropped column" begin
        # Tables.jl would otherwise take the schema from the first row and silently drop every
        # later row's σ_ columns.
        computed = keys(first(plain).val)

        # Both orderings, because the filling works off the first row's names.
        for mixed in ([plain[1], uncertain[1]], [uncertain[1], plain[1]])
            df = DataFrame(mixed)

            @test ncol(df) == 2 * length(computed)
            @test eltype(df.σ_pHtot) === Union{Nothing, Float64}
            @test df.pHtot == [r.pHtot for r in mixed]
            @test count(isnothing, df.σ_pHtot) == 1
        end

        # Column order is values then uncertainties, whichever row came first.
        df = DataFrame([plain[1], uncertain[1]])
        @test names(df)[1:length(computed)] == string.(collect(computed))
        @test all(startswith("σ_"), names(df)[(length(computed) + 1):end])
    end

    @testset "results of different scopes union their columns" begin
        # Same mechanism: a whole-system result carries values a carbon-only one does not.
        carbon = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 10.0)
        whole = whole_system(TA = 2300.0, DIC = 2000.0, temp_c = 10.0, δBT = 39.61)
        df = DataFrame([carbon, whole])

        @test ncol(df) == length(keys(whole.val))
        @test isnothing(df.δBT[1]) && df.δBT[2] == whole.δBT
        @test df.pHtot == [carbon.pHtot, whole.pHtot]
    end

    @testset "edge cases" begin
        # An empty broadcast still knows its columns, because the element type carries them.
        @test size(DataFrame(typeof(first(plain))[])) == (0, length(keys(first(plain).val)))

        # A single result is not a one-row table. `DataFrame([result])` is the spelling.
        @test_throws ArgumentError DataFrame(first(plain))
        @test nrow(DataFrame([first(plain)])) == 1
    end

end
