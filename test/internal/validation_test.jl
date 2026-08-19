# Input validation.
#
# Every case below is input that would otherwise succeed quietly, or fail while naming an
# internal variable rather than the mistake. Grouped by what the user actually did wrong,
# because that is what the error message has to identify.

using Test
using CarbonateCalculator

@testset "Input validation" begin

    @testset "unrecognised keywords" begin
        # Julia routes unmatched keywords into the `kwargs...` sink, so unchecked, `tempc=2.0`
        # is dropped and the calculation runs at the default 25 °C, returning a plausible
        # wrong number.
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, tempc=2.0)
        @test_throws ArgumentError whole_system(TA=2300.0, DIC=2000.0, presbar=100.0)
        @test_throws ArgumentError boron_system(pHtot=8.1, BT=415.7, δBT=39.61, saal=35.0)

        # The message should name the intended argument, not just report a rejection.
        message = try
            carbon_system(TA=2300.0, DIC=2000.0, tempc=2.0)
        catch err
            sprint(showerror, err)
        end
        @test occursin("tempc", message)
        @test occursin("temp_c", message)

        # Retired condition names keep their more specific message.
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, T_in=2.0)
        retired = try
            carbon_system(TA=2300.0, DIC=2000.0, T_out=2.0)
        catch err
            sprint(showerror, err)
        end
        @test occursin("recalculate_at_target_conditions", retired)
    end

    @testset "under-determined" begin
        # Used to surface as `UndefVarError: H not defined in local scope`.
        @test_throws ArgumentError carbon_system()
        @test_throws ArgumentError carbon_system(TA=2300.0)
        @test_throws ArgumentError carbon_system(temp_c=10.0, sal=35.0)

        message = try
            carbon_system(TA=2300.0)
        catch err
            sprint(showerror, err)
        end
        @test !occursin("UndefVarError", message)
        @test occursin("two parameters", message)
    end

    @testset "over-determined" begin
        # Used to silently discard DIC and recompute it: TA=2300, DIC=2000, pHtot=7.0
        # returned DIC = 2445.76 with no indication the input had been dropped.
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, pHtot=7.0)
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, pHtot=8.1)
        @test_throws ArgumentError whole_system(TA=2300.0, DIC=2000.0, pHtot=8.0)
    end

    @testset "one degree of freedom named twice" begin
        # These are the same quantity in different notation, so all but one would have been
        # ignored rather than reconciled.
        @test_throws ArgumentError carbon_system(TA=2300.0, pHtot=8.1, pHsws=8.0)
        @test_throws ArgumentError carbon_system(TA=2300.0, ΩA=3.0, ΩC=5.0)
        @test_throws ArgumentError carbon_system(TA=2300.0, pCO₂=400.0, fCO₂=399.0)
    end

    @testset "boron and isotope routes stay solvable" begin
        # pH is not only a carbon observable: BT with one of BOH₃/BOH₄ fixes it, as does
        # δBT with one of δBOH₃/δBOH₄. The determinacy check must not reject these.
        @test whole_system(BOH₄=100.0, DIC=2000.0).pHtot > 0
        @test whole_system(BOH₃=315.0, DIC=2000.0).pHtot > 0
        @test whole_system(δBOH₄=18.6, DIC=2000.0).pHtot > 0

        # BT and δBT are totals, not constraints, so adding them is not over-determination.
        @test whole_system(BT=415.7, BOH₄=100.0, DIC=2000.0).pHtot > 0
        @test whole_system(δBT=39.61, δBOH₄=18.6, DIC=2000.0).pHtot > 0
    end

    @testset "conditions" begin
        # Unguarded, negative salinity is a DomainError from a sqrt several frames down.
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, sal=-5.0)
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, pres_bar=-10.0)
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, temp_c=NaN)

        message = try
            carbon_system(TA=2300.0, DIC=2000.0, sal=-5.0)
        catch err
            sprint(showerror, err)
        end
        @test occursin("sal", message)

        # Out of range but physical: warn, still compute. Refusing would be worse than
        # extrapolating, and the package should work for a brine.
        @test carbon_system(TA=2300.0, DIC=2000.0, sal=60.0).pHtot > 0
        @test carbon_system(TA=2300.0, DIC=2000.0, temp_c=60.0).pHtot > 0
    end

    @testset "uncertainty names" begin
        # Used to die inside propagate_errors as `type NamedTuple has no field tempc`.
        # Now the name is checked when the solver is built, before any row is computed.
        sigma(name) = CarbonateSystem(:carbon; varying=(:TA, :DIC), varying_errors=(name,))
        @test_throws ArgumentError sigma(:tempc)
        @test_throws ArgumentError sigma(:nonsense)

        # ...while a correctly named one still propagates.
        @test sigma(:TA)(2300.0, 2000.0, 2.0).err.pHtot > 0

        # `errors` is not a parameter of the keyword API: propagation is reached by building
        # a solver with `varying_errors`, so a NamedTuple of σ passed here is a plain typo
        # rather than a silently-ignored request for uncertainties.
        @test_throws ArgumentError carbon_system(TA=2300.0, DIC=2000.0, errors=(TA=2.0,))
    end

    @testset "valid input is unaffected" begin
        # The guards must not have narrowed what the package accepts. Each of these is a
        # legitimate two-parameter call.
        @test carbon_system(TA=2300.0, DIC=2000.0).pHtot > 0
        @test carbon_system(TA=2300.0, pHtot=8.1).DIC > 0
        @test carbon_system(DIC=2000.0, pCO₂=400.0).pHtot > 0
        @test carbon_system(TA=2300.0, ΩC=5.0).pHtot > 0
        @test carbon_system(pHtot=8.1, HCO₃=1800.0).DIC > 0
        @test carbon_system(TA=2300.0, DIC=2000.0, sal=35.0, temp_c=25.0,
                            pres_bar=0.0).pHtot > 0
    end

end
