# Buffer factors computed on demand: ∂X/∂Y at constant Z.
#
# The two hand-written factors in carbon.jl are the reference. They are called by the pipeline
# *before* the concentrations are rescaled, so they work in mol/kg — which is why the reference
# result below is built with `unit = "mol"` rather than passing values off a µmol/kg result.
# Doing the latter returns a plausible number about 1% wrong.

using Test
using CarbonateCalculator

@testset "Gradients and buffer factors" begin

    kw = (TA = 2300.0, DIC = 2000.0, temp_c = 20.0, sal = 35.0, PT = 0.5, SiT = 3.0)
    result = carbon_system(; kw...)

    @testset "agrees with the hand-written factors" begin
        in_mol = carbon_system(; TA = 2300e-6, DIC = 2000e-6, temp_c = 20.0, sal = 35.0,
                               PT = 0.5e-6, SiT = 3e-6, unit = "mol")
        v = in_mol.val

        # Reached through `Carbon` rather than the package: both are kept for reference but no
        # longer exported, so they do not reach `CarbonateCalculator` via `using .Carbon`.
        reference = CarbonateCalculator.Carbon.calc_revelle_factor(
            v.TA, v.DIC, v.BT, v.PT, v.SiT, v.ST, v.FT, v.H2ST, v.NH4T, in_mol.Ks)
        @test revelle_factor(in_mol) ≈ reference rtol=1e-12

        capacity = CarbonateCalculator.Carbon.calc_buffer_capacity(
            v.pHtot, v.DIC, v.BT, v.PT, v.SiT, v.ST, v.FT, v.H2ST, v.NH4T, in_mol.Ks)
        @test calc_gradient(in_mol, :TA, :pHtot; constant = :DIC) ≈ capacity rtol=1e-12
    end

    @testset "the held constant is half the definition" begin
        at_TA = calc_gradient(result, :pCO₂, :DIC; constant = :TA)
        at_pH = calc_gradient(result, :pCO₂, :DIC; constant = :pHtot)

        # Not a small difference: roughly a factor of ten. Naming numerator and denominator
        # alone would not have specified which of these was meant.
        @test at_TA > 5 * at_pH
        @test at_TA > 0 && at_pH > 0
    end

    @testset "relative is the gradient normalised" begin
        gradient = calc_gradient(result, :pCO₂, :DIC; constant = :TA)
        relative = calc_relative_gradient(result, :pCO₂, :DIC; constant = :TA)
        @test relative ≈ gradient * result.DIC / result.pCO₂ rtol=1e-14

        # fCO₂ and pCO₂ differ by a factor fixed at constant T, S and P, so it cancels.
        @test revelle_factor(result) ≈
              calc_relative_gradient(result, :pCO₂, :DIC; constant = :TA) rtol=1e-14

        # Dimensionless, so the reporting unit cannot matter.
        in_mol = carbon_system(; TA = 2300e-6, DIC = 2000e-6, temp_c = 20.0, sal = 35.0,
                               PT = 0.5e-6, SiT = 3e-6, unit = "mol")
        @test revelle_factor(result) ≈ revelle_factor(in_mol) rtol=1e-10
    end

    @testset "ill-posed requests are refused" begin
        # Varying and holding the same quantity.
        @test_throws "both varied and held" calc_gradient(result, :pCO₂, :DIC; constant = :DIC)

        # pCO₂ and fCO₂ are one quantity expressed two ways, so they cannot be an independent
        # pair. This comes from the solver's own determinacy check rather than a second list.
        @test_throws ArgumentError calc_gradient(result, :pHtot, :fCO₂; constant = :pCO₂)

        # A name that is not a computed value at all.
        @test_throws Exception calc_gradient(result, :pCO₂, :DIC; constant = :not_a_parameter)
    end

    @testset "works for every kind of result" begin
        # Boron in scope: the constraint set comes from PARAMETER_GROUPS, so the boron
        # parameters are stripped too and cannot over-determine the re-solve.
        whole = whole_system(; kw..., δBT = 39.61)
        @test revelle_factor(whole) > 0
        @test calc_gradient(whole, :pHtot, :DIC; constant = :TA) < 0   # more DIC, lower pH

        # A result already carried to its collection conditions.
        collected = at_collection_conditions(result; temp_c = 2.0)
        @test revelle_factor(collected) != revelle_factor(result)

        # A result whose measured pair was not TA/DIC, so the held value is one it derived.
        from_pH = carbon_system(; pHtot = 8.1, TA = 2300.0, temp_c = 20.0, sal = 35.0)
        @test revelle_factor(from_pH) > 0
    end

    @testset "uncertainty on a gradient" begin
        measured = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC), TA = 2300.0,
                                   DIC = 2000.0, temp_c = 20.0, sal = 35.0,
                                   PT = 0.5, SiT = 3.0)(2.0, 2.0)
        σ = calc_gradient_uncertainty(measured, :pCO₂, :DIC; constant = :TA)

        # The guard that matters. Freezing the evaluation point and the held constant at their
        # solved values returns exactly zero here, because the carbonate constraints are
        # stripped before the re-solve and the perturbation never lands. Nonzero is the whole
        # correctness claim.
        @test σ > 0

        # Validated against finite differences, because nested AD fails silently when it fails.
        # `calc_gradient` off a TA/DIC result is a function of those two inputs, so a central
        # difference in each gives the same sensitivities the second-order pass computes.
        gradient_at(TA, DIC) = calc_gradient(
            carbon_system(; TA = TA, DIC = DIC, temp_c = 20.0, sal = 35.0, PT = 0.5, SiT = 3.0),
            :pCO₂, :DIC; constant = :TA)

        hTA, hDIC = 0.23, 0.20
        dTA = (gradient_at(2300.0 + hTA, 2000.0) - gradient_at(2300.0 - hTA, 2000.0)) / 2hTA
        dDIC = (gradient_at(2300.0, 2000.0 + hDIC) - gradient_at(2300.0, 2000.0 - hDIC)) / 2hDIC
        @test σ ≈ sqrt((dTA * 2.0)^2 + (dDIC * 2.0)^2) rtol=1e-5

        # Bigger input uncertainties can only widen it, and it scales linearly at first order.
        wider = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC), TA = 2300.0,
                                DIC = 2000.0, temp_c = 20.0, sal = 35.0,
                                PT = 0.5, SiT = 3.0)(4.0, 4.0)
        @test calc_gradient_uncertainty(wider, :pCO₂, :DIC; constant = :TA) ≈ 2σ rtol=1e-9

        # A σ that reaches the answer only through the equilibrium constants still propagates.
        via_constants = CarbonateSystem(:carbon; varying_errors = (:temp_c,), TA = 2300.0,
                                        DIC = 2000.0, temp_c = 20.0, sal = 35.0)(0.1)
        @test calc_gradient_uncertainty(via_constants, :pCO₂, :DIC; constant = :TA) > 0

        # Nothing to propagate.
        @test isnothing(calc_gradient_uncertainty(result, :pCO₂, :DIC; constant = :TA))

        # Ill-posed requests are refused on this path too, not just the value one.
        @test_throws "both varied and held" calc_gradient_uncertainty(
            measured, :pCO₂, :DIC; constant = :DIC)
    end

    @testset "uncertainty at collection conditions runs through both stages" begin
        # A non-conservative measured pair, so stage 1 does real work and the chain matters:
        # `input_errors` names pHtot, which the collection state does not carry at all.
        base = (sal = 35.0, PT = 0.5, SiT = 3.0)
        measured = CarbonateSystem(:carbon; varying_errors = (:pHtot, :TA), pHtot = 8.1,
                                   TA = 2300.0, temp_c = 20.0, base...)(0.005, 2.0)

        # 10 °C deliberately, not 2 °C: `which_K` switches parameterisation at exactly 2 °C
        # for this salinity, and a derivative across that step is meaningless.
        collected = at_collection_conditions(measured; temp_c = 10.0, σ_temp_c = 0.5)

        @test !isnothing(collected.source)
        @test collected.source.condition_errors == (target_temp_c = 0.5,)
        @test isnothing(collected.source.inputs.pHtot) == false   # the measurement kept its pH

        σ = calc_gradient_uncertainty(collected, :pCO₂, :DIC; constant = :TA)
        @test σ > 0

        # Validated against finite differences over the whole two-stage chain, perturbing what
        # was actually measured plus the collection temperature.
        chain(pH, TA, t) = calc_gradient(
            at_collection_conditions(carbon_system(; pHtot = pH, TA = TA, temp_c = 20.0,
                                                   base...); temp_c = t),
            :pCO₂, :DIC; constant = :TA)

        hpH, hTA, ht = 1e-5, 0.23, 2e-3
        dpH = (chain(8.1 + hpH, 2300.0, 10.0) - chain(8.1 - hpH, 2300.0, 10.0)) / 2hpH
        dTA = (chain(8.1, 2300.0 + hTA, 10.0) - chain(8.1, 2300.0 - hTA, 10.0)) / 2hTA
        dt = (chain(8.1, 2300.0, 10.0 + ht) - chain(8.1, 2300.0, 10.0 - ht)) / 2ht
        @test σ ≈ sqrt((dpH * 0.005)^2 + (dTA * 2.0)^2 + (dt * 0.5)^2) rtol=1e-5

        # The collection temperature's own uncertainty is a real contributor, not a rounding
        # term — dropping it would understate σ here by about a third.
        without_condition_σ = at_collection_conditions(measured; temp_c = 10.0)
        @test calc_gradient_uncertainty(without_condition_σ, :pCO₂, :DIC; constant = :TA) < σ

        # A collection state with nothing uncertain anywhere has no σ to report.
        plain = at_collection_conditions(
            carbon_system(; TA = 2300.0, DIC = 2000.0, temp_c = 20.0, base...); temp_c = 10.0)
        @test isnothing(calc_gradient_uncertainty(plain, :pCO₂, :DIC; constant = :TA))
    end

    @testset "a relative gradient's uncertainty is not the scaled one" begin
        measured = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC), TA = 2300.0,
                                   DIC = 2000.0, temp_c = 20.0, sal = 35.0,
                                   PT = 0.5, SiT = 3.0)(2.0, 2.0)

        relative = calc_gradient_uncertainty(measured, :fCO₂, :DIC; constant = :TA,
                                             relative = true)
        raw = calc_gradient_uncertainty(measured, :fCO₂, :DIC; constant = :TA)

        # `(∂X/∂Y)·(Y/X)` normalises by two quantities that are themselves uncertain, so
        # scaling the raw σ afterwards misses their contribution. Materially different, not a
        # rounding difference.
        naive = raw * measured.DIC / measured.fCO₂
        @test relative > 0
        @test !isapprox(relative, naive; rtol = 0.05)
    end

    @testset "with_gradient adds a value without disturbing the result" begin
        added = with_gradient(result, :revelle, :fCO₂, :DIC; constant = :TA, relative = true)

        @test added.revelle ≈ revelle_factor(result) rtol=1e-14
        @test keys(added.val) == (keys(result.val)..., :revelle)
        @test !haskey(result.val, :revelle)          # the original is untouched
        @test added.pHtot == result.pHtot            # everything else survives
        @test added.system == result.system
        @test added.source === result.source

        # Absolute as well as relative.
        raw = with_gradient(result, :dpCO2_dDIC, :pCO₂, :DIC; constant = :TA)
        @test raw.dpCO2_dDIC ≈ calc_gradient(result, :pCO₂, :DIC; constant = :TA) rtol=1e-14

        # A name already in use would silently shadow a computed value.
        @test_throws "already has a value called" with_gradient(
            result, :pHtot, :pCO₂, :DIC; constant = :TA)

        # The source survives, so a carried result keeps refusing a second carry.
        collected = with_gradient(at_collection_conditions(result; temp_c = 2.0),
                                  :revelle, :fCO₂, :DIC; constant = :TA, relative = true)
        @test !isnothing(collected.source)
        @test_throws ArgumentError at_collection_conditions(collected; temp_c = 4.0)
    end

    @testset "with_gradient carries the uncertainty across too" begin
        measured = CarbonateSystem(:carbon; varying_errors = (:TA, :DIC), TA = 2300.0,
                                   DIC = 2000.0, temp_c = 20.0, sal = 35.0,
                                   PT = 0.5, SiT = 3.0)(2.0, 2.0)
        added = with_gradient(measured, :revelle, :fCO₂, :DIC; constant = :TA,
                              relative = true)

        @test added.err.revelle ≈ calc_gradient_uncertainty(measured, :fCO₂, :DIC;
                                                            constant = :TA,
                                                            relative = true) rtol=1e-14
        # `val` and `err` stay key-for-key aligned, which is what the σ_ columns assume.
        @test keys(added.val) == keys(added.err)

        # No measurement uncertainties: the value lands, `err` stays absent.
        plain = with_gradient(result, :revelle, :fCO₂, :DIC; constant = :TA, relative = true)
        @test plain.revelle > 0
        @test isnothing(plain.err)

        # A collection state gets both, through the two-stage chain.
        carried = with_gradient(at_collection_conditions(measured; temp_c = 10.0,
                                                         σ_temp_c = 0.5),
                                :revelle, :fCO₂, :DIC; constant = :TA, relative = true)
        @test carried.revelle > 0
        @test carried.err.revelle > 0
        @test keys(carried.val) == keys(carried.err)
    end

end
