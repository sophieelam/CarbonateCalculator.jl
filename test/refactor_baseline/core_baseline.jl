# Re-capturable bitwise snapshot of the current API.
#
#   julia --project=test test/refactor_baseline/core_baseline.jl capture
#   julia --project=test test/refactor_baseline/core_baseline.jl check
#
# Separate from the two files beside it, which cannot do this job:
#   * `capture_baseline.jl` is frozen — every case uses the retired `T_in`/`T_out`
#     arguments, so it now raises and is kept only as the record of how
#     `baseline_two_condition.tsv` was produced.
#   * `check_against_baseline.jl` compares against that frozen TSV with known-bug
#     allowances. Still valuable, still run, but its 18 cases are all condition changes
#     and would not notice a units or method-matrix regression.
#
# Stage C needs the opposite: a snapshot that can be *re-taken* before each step, over a
# wide sweep, compared **bitwise**. Bitwise on purpose — restructuring changes the order of
# operations, and a last-ulp shift is exactly what we want surfaced and explained rather
# than absorbed by a tolerance.

using CarbonateCalculator
using Printf

const OUT = joinpath(@__DIR__, "core_baseline.tsv")
const CONDITIONS = ((25.0, 35.0, 0.0), (2.0, 35.0, 400.0), (10.0, 30.0, 100.0))

"Every case, as `name => kwargs`. Deterministic, and wider than the test suite."
function cases()
    out = Pair{String,Any}[]

    pairs = ((:TA => 2300.0, :DIC => 2000.0), (:TA => 2300.0, :pHtot => 8.1),
             (:DIC => 2000.0, :pHtot => 8.1), (:TA => 2300.0, :pCO₂ => 400.0),
             (:DIC => 2000.0, :pCO₂ => 400.0), (:TA => 2300.0, :fCO₂ => 400.0),
             (:TA => 2300.0, :CO₃ => 200.0), (:TA => 2300.0, :HCO₃ => 1800.0),
             (:DIC => 2000.0, :CO₃ => 200.0), (:pHtot => 8.1, :CO₃ => 200.0),
             (:TA => 2300.0, :ΩC => 5.0), (:TA => 2300.0, :ΩA => 3.0))
    for (temp_c, sal, pres_bar) in CONDITIONS, (a, b) in pairs
        push!(out, "pair_$(a[1])_$(b[1])_T$(temp_c)_P$(pres_bar)" =>
                   (; a[1] => a[2], b[1] => b[2], temp_c, sal, pres_bar))
    end

    for K_method in ("default", "Roy 1993", "GP 1989", "Hansson 1973", "DM 1987",
                     "HM 1973", "Mehrbach 1973 A", "Mehrbach 1973 B", "Millero 1979",
                     "CW 2003", "Lueker 2000", "MPM 2002", "Millero 2002", "Millero 2006",
                     "Millero 2010", "Waters 2014", "SB 2020", "Papadimitriou 2018",
                     "Sulpis 2020", "MyAMI"),
        (temp_c, sal, pres_bar) in CONDITIONS
        push!(out, "K_$(K_method)_T$(temp_c)_P$(pres_bar)" =>
                   (TA = 2300.0, DIC = 2000.0, temp_c, sal, pres_bar, K_method))
    end

    for KSO4_method in ("Dickson", "Khoo", "WM13"), BT_method in ("Uppstrom", "Lee", "KSK18")
        push!(out, "KSO4_$(KSO4_method)_BT_$(BT_method)" =>
                   (TA = 2300.0, DIC = 2000.0, KSO4_method, BT_method))
    end
    for KF_method in ("Dickson", "Perez"), KNH3_method in ("Millero", "Clegg"),
        Ca_method in ("default", "Culkin", "RT67")
        push!(out, "KF_$(KF_method)_KNH3_$(KNH3_method)_Ca_$(Ca_method)" =>
                   (TA = 2300.0, DIC = 2000.0, KF_method, KNH3_method, Ca_method))
    end

    # Units matter here: the unit rescaling is one of the blocks being consolidated.
    for (unit, TA, DIC) in (("mol", 2.3e-3, 2.0e-3), ("mmol", 2.3, 2.0),
                            ("umol", 2300.0, 2000.0), ("nmol", 2.3e6, 2.0e6))
        push!(out, "unit_$unit" => (; TA, DIC, unit, PT = 1.0, SiT = 10.0,
                                    H2ST = 3.0, NH4T = 2.0))
    end
    push!(out, "nutrients" => (TA = 2300.0, DIC = 2000.0, PT = 2.5, SiT = 50.0,
                               H2ST = 5.0, NH4T = 3.0))
    push!(out, "paleo_CaMg" => (TA = 2300.0, DIC = 2000.0, Ca = 0.0206, Mg = 0.0323))
    push!(out, "explicit_totals" => (TA = 2300.0, DIC = 2000.0, ST = 0.028,
                                     FT = 7.0e-5, BT = 420.0))
    return out
end

"Boron arguments added when a case is also run through whole_system."
const BORON = (δBT = 39.61, alphaB = 1.0272)

"""
Cases for the two boron-only entry points, which take a different argument set.

`boron_system` and `boron_isotopes` are covered because phase 5b folds them onto the same
stages as the cores, and without a snapshot that refactor would have no safety net.
"""
function boron_cases()
    out = Pair{String,Any}[]
    for (temp_c, sal, pres_bar) in CONDITIONS
        for BT in (415.7, 300.0), pHtot in (8.1, 7.6)
            out = push!(out, "b_pH$(pHtot)_BT$(BT)_T$(temp_c)" =>
                             (; pHtot, BT, δBT = 39.61, temp_c, sal, pres_bar))
        end
        push!(out, "b_BOH4_T$(temp_c)" =>
                   (BT = 415.7, BOH₄ = 100.0, δBT = 39.61, temp_c, sal, pres_bar))
        push!(out, "b_dBOH4_T$(temp_c)" =>
                   (BT = 415.7, δBOH₄ = 18.6, δBT = 39.61, temp_c, sal, pres_bar))
        push!(out, "b_alphaB_T$(temp_c)" =>
                   (pHtot = 8.1, BT = 415.7, δBT = 39.61, alphaB = 1.0272,
                    temp_c, sal, pres_bar))
    end
    for K_method in ("default", "Lueker 2000", "MyAMI")
        push!(out, "b_K_$(K_method)" =>
                   (pHtot = 8.1, BT = 415.7, δBT = 39.61, K_method = K_method))
    end
    return out
end

"Every finite numeric output of both systems for one case, as (system, key, value)."
function evaluate(kwargs)
    rows = Tuple{String,String,Float64}[]
    for (system, fn, extra) in (("carbon", carbon_system, NamedTuple()),
                                ("whole", whole_system, BORON))
        result = try
            fn(; kwargs..., extra...)
        catch
            push!(rows, (system, "THREW", NaN))
            continue
        end
        for (key, value) in pairs(result.val)
            value isa Real && isfinite(value) &&
                push!(rows, (system, string(key), Float64(value)))
        end
    end
    return rows
end

"""
Same, for the boron-only entry points.

These return a bare NamedTuple rather than a `CarbonateResult`, so there is no `.val`.
"""
function evaluate_boron(kwargs)
    rows = Tuple{String,String,Float64}[]
    for (system, fn) in (("boron_system", boron_system),
                         ("boron_isotopes", boron_isotopes))
        result = try
            fn(; kwargs...)
        catch
            push!(rows, (system, "THREW", NaN))
            continue
        end
        for (key, value) in pairs(result)
            value isa Real && isfinite(value) &&
                push!(rows, (system, string(key), Float64(value)))
        end
    end
    return rows
end

"Every (case, system, key, value) row the snapshot covers."
function all_rows()
    rows = Tuple{String,String,String,Float64}[]
    for (name, kwargs) in cases(), (system, key, value) in evaluate(kwargs)
        push!(rows, (name, system, key, value))
    end
    for (name, kwargs) in boron_cases(), (system, key, value) in evaluate_boron(kwargs)
        push!(rows, (name, system, key, value))
    end
    return rows
end

function capture(path)
    rows = all_rows()
    open(path, "w") do io
        println(io, "case\tsystem\tkey\tvalue")
        for (name, system, key, value) in rows
            # %.17g round-trips a Float64 exactly, so the file is a bitwise record.
            @printf(io, "%s\t%s\t%s\t%.17g\n", name, system, key, value)
        end
    end
    println("captured $(length(rows)) values from " *
            "$(length(cases()) + length(boron_cases())) cases to $(basename(path))")
end

function check(path)
    isfile(path) || error("no snapshot at $path — run `capture` first")

    saved = Dict{Tuple{String,String,String},Float64}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        case, system, key, value = split(line, '\t')
        saved[(case, system, key)] = parse(Float64, value)
    end

    matched, differing, appeared = 0, Tuple{String,Float64}[], 0
    for (name, system, key, value) in all_rows()
        entry = (name, system, key)
        if !haskey(saved, entry)
            appeared += 1
            continue
        end
        before = pop!(saved, entry)
        if before === value || (isnan(before) && isnan(value))
            matched += 1
        else
            rel = before == 0 ? abs(value) : abs((value - before) / before)
            push!(differing, ("$name / $system / $key: $before -> $value", rel))
        end
    end

    println("matched bitwise : $matched")
    println("differing       : $(length(differing))")
    println("absent now      : $(length(saved))")
    println("new now         : $appeared")

    if !isempty(differing)
        sort!(differing, by = last, rev = true)
        println("\nlargest relative differences:")
        for (text, rel) in first(differing, min(25, length(differing)))
            @printf("  %-11.3e  %s\n", rel, text)
        end
    end
    for entry in first(sort(collect(keys(saved))), min(10, length(saved)))
        println("  absent now: ", entry)
    end

    return isempty(differing) && isempty(saved) && appeared == 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "check" : ARGS[1]
    mode == "capture" ? capture(OUT) :
    mode == "check"   ? exit(check(OUT) ? 0 : 1) :
    error("usage: core_baseline.jl [capture|check]")
end
