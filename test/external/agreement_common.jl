# Shared harness for the external-calculator agreement tests (compare_*.jl).
#
# Each test hands a baseline CSV (CO2SYS-MATLAB column names; see
# generate_pyco2sys_baseline.py and compare_MATLABv3_2_0.m) and one of the package's entry
# points to `run_agreement`, which solves every row, compares everything the baseline holds
# and the package computes, and fails on any difference beyond fit-for-purpose tolerances.
#
# Principles, so that this test can fail when the chemistry is wrong:
#
# * No pattern-based exemptions. A row that disagrees is a failure unless its K1/K2 option
#   is listed in DOCUMENTED_DIFFERENCES with the reason, which is reported separately (and
#   should match the README).
# * Tolerances are CarbCalcComp's (3C) fit-for-purpose ones, in the spirit of the GOA-ON
#   "climate" goals: pH 3e-4, TA and DIC 0.2 µmol/kg, pCO2 0.05 %, CO3 and Ω 0.1 %. PyCO2SYS
#   v1.8.3 and CO2SYS-MATLAB v3.2 agree with each other to 1.4e-7 in pH on every row of
#   these baselines, so anything near these tolerances is a real difference.
# * Values at the output conditions and the equilibrium constants are compared too: the
#   pressure correction is only visible there, or at a pressurised input.
#
# Every failure is written to results/failures_<name>.csv for inspection.
module Agreement

using Test, CSV, DataFrames, Printf
using CarbonateCalculator

export run_agreement

# CO2SYS's K1K2CONSTANTS options as K_settings. Options 6-8 bring the companion constants
# CO2SYS switches to with them (GEOSECS, Peng et al., freshwater), and those win over the
# BORON option, as in CO2SYS.
const K1K2 = Dict(
    1 => (K1_K2_method = "Roy 1993",),
    2 => (K1_K2_method = "GP 1989",),
    3 => (K1_K2_method = "Hansson 1973",),
    4 => (K1_K2_method = "DM 1987",),
    5 => (K1_K2_method = "HM 1973",),
    6 => (K1_K2_method = "Mehrbach 1973", BT_method = "Culkin", Ca_method = "Culkin"),
    7 => (K1_K2_method = "Mehrbach 1973", fH_method = "Peng 1987", BT_method = "Culkin",
          Ca_method = "Culkin"),
    8 => (K1_K2_method = "Millero 1979", fH_method = "freshwater", BT_method = "freshwater",
          pressure_method = "Millero 1983"),
    9 => (K1_K2_method = "CW 2003",),
    10 => (K1_K2_method = "Lueker 2000",),
    11 => (K1_K2_method = "MPM 2002",),
    12 => (K1_K2_method = "Millero 2002",),
    13 => (K1_K2_method = "Millero 2006",),
    14 => (K1_K2_method = "Millero 2010",),
    15 => (K1_K2_method = "Waters 2014",),
    16 => (K1_K2_method = "Sulpis 2020",),
    17 => (K1_K2_method = "SB 2020",),
)
const KSO4 = Dict(1 => "Dickson", 2 => "Khoo", 3 => "WM13")
const BORON = Dict(1 => "Uppstrom", 2 => "Lee")
const KF = Dict(1 => "Dickson", 2 => "Perez")
# Input-type codes as in PyCO2SYS (6 = CO3, 7 = HCO3); 3 = pH on the row's input scale.
const PAR = Dict(1 => :TA, 2 => :DIC, 3 => :pH, 4 => :pCO₂, 5 => :fCO₂, 6 => :CO₃, 7 => :HCO₃,
                 8 => :CO₂)
const PH_FIELD = Dict(1 => :pHtot, 2 => :pHsws, 3 => :pHfree, 4 => :pHNBS)

# Known, intended differences, by K1/K2 option. Keep in step with the README.
# (README, "CO2SYS changes more than this under options 6-8".)
const DOCUMENTED_DIFFERENCES = Dict(
    6 => "GEOSECS parts of CO2SYS not implemented (README, options 6-8)",
    7 => "Peng et al. parts of CO2SYS not implemented (README, options 6-8)",
    8 => "freshwater with nutrients: CO2SYS zeroes KB, KP, KSi, KH2S, KNH4 (README, option 8)",
)

"""
Is a disagreement on this row a documented difference? Options 6 and 7 always; option 8
only where phosphate, silicate, sulfide or ammonia is given, because CO2SYS sets their
constants to zero for fresh water (so they add no alkalinity) and this package does not.
"""
documented(k1k2, row) = k1k2 in (6, 7) ||
    (k1k2 == 8 && any(row[c] > 0 for c in (:PO4, :SI, :H2S, :NH4)))

# (label, result field or function, baseline column, tolerance). Tolerances are absolute
# (atol) or relative (rtol), as in 3C's config/variables.yaml.
const AT_INPUT = [
    ("TA", :TA, "TAlk", (atol = 0.2,)),
    ("DIC", :DIC, "TCO2", (atol = 0.2,)),
    ("pH (input scale)", :pH, "pHin", (atol = 3e-4,)),
    ("pH total", :pHtot, "pHinTOTAL", (atol = 3e-4,)),
    ("pCO2", :pCO₂, "pCO2in", (rtol = 5e-4,)),
    ("fCO2", :fCO₂, "fCO2in", (rtol = 5e-4,)),
    ("CO2", :CO₂, "CO2in", (rtol = 5e-4,)),
    ("HCO3", :HCO₃, "HCO3in", (rtol = 1e-3,)),
    ("CO3", :CO₃, "CO3in", (rtol = 1e-3,)),
    ("Ω calcite", :ΩC, "OmegaCAin", (rtol = 1e-3,)),
    ("Ω aragonite", :ΩA, "OmegaARin", (rtol = 1e-3,)),
    ("Revelle factor", :revelle, "RFin", (rtol = 1e-2,)),
]
const AT_OUTPUT = [
    ("pH out (input scale)", :pH, "pHout", (atol = 3e-4,)),
    ("pH total out", :pHtot, "pHoutTOTAL", (atol = 3e-4,)),
    ("fCO2 out", :fCO₂, "fCO2out", (rtol = 5e-4,)),
    ("CO2 out", :CO₂, "CO2out", (rtol = 5e-4,)),
    ("HCO3 out", :HCO₃, "HCO3out", (rtol = 1e-3,)),
    ("CO3 out", :CO₃, "CO3out", (rtol = 1e-3,)),
    ("Ω calcite out", :ΩC, "OmegaCAout", (rtol = 1e-3,)),
    ("Ω aragonite out", :ΩA, "OmegaARout", (rtol = 1e-3,)),
]
# Constants at the input conditions; compared only for total-scale rows, because the
# package reports them on the total scale and CO2SYS on the input scale (KS, KF free).
const CONSTANTS = [
    ("K0", :K0, "K0input", (rtol = 5e-4,)),
    ("K1", :K1, "K1input", (rtol = 1e-4,)),
    ("K2", :K2, "K2input", (rtol = 1e-4,)),
    ("KB", :KB, "KBinput", (rtol = 1e-4,)),
    ("KW", :KW, "KWinput", (rtol = 1e-4,)),
    ("KS", :KS, "KSinput", (rtol = 1e-3,)),
    ("KF", :KF, "KFinput", (rtol = 1e-3,)),
    ("KP1", :KP1, "KP1input", (rtol = 1e-4,)),
    ("KP2", :KP2, "KP2input", (rtol = 1e-4,)),
    ("KP3", :KP3, "KP3input", (rtol = 1e-4,)),
    ("KSi", :KSi, "KSiinput", (rtol = 1e-4,)),
    ("KH2S", :KH2S, "KH2Sinput", (rtol = 1e-4,)),
    ("KNH3", :KNH3, "KNH3input", (rtol = 1e-4,)),
    ("Ksp calcite", :KspC, "KSPCinput", (rtol = 1e-3,)),
    ("Ksp aragonite", :KspA, "KSPAinput", (rtol = 1e-3,)),
]

unpack(x) = Float64(x isa Tuple ? first(x) : x)

"The package's value, or NaN if this result doesn't provide it (then counted as a failure)."
function value(res, field, scale)
    try
        field === :pH && return unpack(getproperty(res, PH_FIELD[scale]))
        field === :revelle && return Float64(CarbonateCalculator.revelle_factor(res))
        return unpack(getproperty(res, field))
    catch
        return NaN
    end
end

"Error measure and verdicts: (error, within tolerance, identical to ~1e-6)."
function judge(calc, ref, tol)
    if haskey(tol, :atol)
        err = abs(calc - ref)
        return err, err <= tol.atol, err <= 1e-6
    end
    err = abs(calc - ref) / max(abs(ref), 1e-300)
    return err, err <= tol.rtol, err <= 1e-6
end

"""
    run_agreement(name, path, solve)

Compare `solve(; inputs..., env...)` (e.g. `whole_system`) with every row of the baseline at
`path`. Returns the failures as a DataFrame and runs `@test` on their number.
"""
function run_agreement(name, path, solve)
    df = CSV.read(path, DataFrame)
    cols = Set(names(df))
    failures = DataFrame(row = Int[], K1K2 = Int[], pHscale = Int[], pair = String[],
                         temp_c = Float64[], sal = Float64[], pres_dbar = Float64[],
                         variable = String[], package = Float64[], reference = Float64[],
                         error = Float64[], tolerance = String[])
    counts = Dict(:compared => 0, :identical => 0, :agree => 0, :fail => 0,
                  :documented => 0, :errors => 0, :skipped => 0)
    documented_n = Dict{Int, Int}()
    out_unsupported = false

    for (i, row) in enumerate(eachrow(df))
        k1k2 = Int(row.K1K2CONSTANTS)
        # The freshwater constants (option 8) are only meaningful at S = 0, the others at S > 0
        if (k1k2 == 8) != (row.SAL == 0)
            counts[:skipped] += 1
            continue
        end
        scale = Int(row.pHSCALEIN)
        boron = "BORON" in cols ? Int(row.BORON) : 1
        settings = merge((KF_method = KF[Int(row.KFCONSTANT)],
                          KSO4_method = KSO4[Int(row.KSO4CONSTANT)],
                          BT_method = BORON[boron],
                          KNH3_method = "Clegg"),   # Clegg and Whitfield (1995), as CO2SYS v3
                         K1K2[k1k2])
        env = (temp_c = row.TEMPIN, sal = row.SAL, pres_bar = row.PRESIN / 10, PT = row.PO4,
               SiT = row.SI, H2ST = row.H2S, NH4T = row.NH4, unit = "umol",
               K_settings = settings)

        # GEOSECS and Peng have no fugacity correction (pCO2 = fCO2), so a pCO2 input is
        # given as fCO2, as CO2SYS does.
        geosecs = k1k2 in (6, 7)
        inputs = Dict{Symbol, Any}()
        pair = String[]
        for (code, val) in ((Int(row.PAR1TYPE), row.PAR1), (Int(row.PAR2TYPE), row.PAR2))
            sym = PAR[code]
            push!(pair, string(sym))
            if sym === :pH
                inputs[PH_FIELD[scale]] = val
            elseif sym === :pCO₂
                inputs[geosecs ? :fCO₂ : :pCO₂] = row.pCO2in
            elseif sym === :fCO₂
                inputs[:fCO₂] = row.fCO2in
            else
                inputs[sym] = val
            end
        end
        pair_name = join(pair, "+")
        is_documented = documented(k1k2, row)

        function record(label, calc, ref, tol)
            (ismissing(ref) || isnan(ref)) && return
            counts[:compared] += 1
            err, ok, same = isfinite(calc) ? judge(calc, ref, tol) : (Inf, false, false)
            if ok
                counts[:agree] += 1
                same && (counts[:identical] += 1)
            elseif is_documented
                counts[:documented] += 1
                documented_n[k1k2] = get(documented_n, k1k2, 0) + 1
            else
                counts[:fail] += 1
                push!(failures, (i, k1k2, scale, pair_name, row.TEMPIN, row.SAL, row.PRESIN,
                                 label, calc, ref, err, string(tol)))
            end
        end

        res = try
            solve(; inputs..., env...)
        catch e
            counts[:errors] += 1
            push!(failures, (i, k1k2, scale, pair_name, row.TEMPIN, row.SAL, row.PRESIN,
                             "error: " * first(sprint(showerror, e), 120), NaN, NaN, Inf, ""))
            continue
        end

        for (label, field, col, tol) in AT_INPUT
            col in cols && record(label, value(res, field, scale), df[i, col], tol)
        end
        if scale == 1 && hasproperty(res, :Ks)
            for (label, field, col, tol) in CONSTANTS
                (col in cols && hasproperty(res.Ks, field)) || continue
                ref = df[i, col]
                if k1k2 == 8 && !ismissing(ref) && ref == 0
                    # CO2SYS's freshwater zeroing (option 8), whatever the totals
                    counts[:compared] += 1; counts[:documented] += 1
                    documented_n[8] = get(documented_n, 8, 0) + 1
                    continue
                end
                record(label, Float64(getproperty(res.Ks, field)), ref, tol)
            end
        end
        if "TEMPOUT" in cols && !out_unsupported
            res_out = try
                CarbonateCalculator.adjust(res, row.TEMPOUT, row.PRESOUT / 10)
            catch e
                if e isa MethodError
                    out_unsupported = true   # this entry point's result can't be adjusted
                    nothing
                else
                    counts[:errors] += 1
                    push!(failures, (i, k1k2, scale, pair_name, row.TEMPIN, row.SAL,
                                     row.PRESIN, "error (output conditions): " *
                                     first(sprint(showerror, e), 120), NaN, NaN, Inf, ""))
                    nothing
                end
            end
            if res_out !== nothing
                for (label, field, col, tol) in AT_OUTPUT
                    col in cols && record(label, value(res_out, field, scale), df[i, col], tol)
                end
            end
        end
    end

    # --- Report -------------------------------------------------------------------------
    println("\n", "═"^72)
    println("  $name: $(nrow(df)) baseline rows, $(counts[:compared]) values compared")
    println("  ", "─"^68)
    @printf("  within tolerance        %10d   (of which ≤ 1e-6: %d)\n", counts[:agree],
            counts[:identical])
    @printf("  FAILURES                %10d\n", counts[:fail])
    @printf("  errors                  %10d\n", counts[:errors])
    for (k, n) in sort(collect(documented_n))
        @printf("  documented (option %2d)  %10d   %s\n", k, n, DOCUMENTED_DIFFERENCES[k])
    end
    counts[:skipped] > 0 &&
        @printf("  skipped                 %10d   (freshwater option off S = 0, or vice versa)\n",
                counts[:skipped])
    out_unsupported && println("  output conditions not compared: this entry point's ",
                               "result does not support `adjust`")
    if nrow(failures) > 0
        file = joinpath(dirname(path), "failures_$(replace(name, r"[^A-Za-z0-9]+" => "_")).csv")
        CSV.write(file, failures)
        println("  ", "─"^68)
        println("  Largest failures by variable (all in $(basename(file))):")
        worst = combine(groupby(failures, :variable), nrow => :n, :error => maximum => :max_error)
        for r in eachrow(first(sort(worst, :n, rev = true), 15))
            @printf("    %-28s %7d   max error %.3g\n", r.variable, r.n, r.max_error)
        end
        byopt = combine(groupby(failures, :K1K2), nrow => :n)
        println("  By K1/K2 option: ", join(["$(r.K1K2): $(r.n)" for r in eachrow(sort(byopt, :K1K2))], ", "))
    end
    println("═"^72, "\n")

    @test counts[:fail] == 0
    @test counts[:errors] == 0
    return failures
end

end # module Agreement
