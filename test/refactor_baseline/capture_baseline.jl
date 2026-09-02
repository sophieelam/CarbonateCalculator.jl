# Phase 0 - capture current two-condition behaviour before the refactor removes it.
# Writes a TSV of case_id / key / value for every numeric field of `.val`.
# Cases deliberately include the four known bugs so that post-refactor differences
# can be classified as "expected bug fix" vs "unexplained regression".
#
#   julia --project=. test/refactor_baseline/capture_baseline.jl [outfile]
#
# `baseline_two_condition.tsv` is the pre-refactor reference and must not be regenerated -
# it cannot be reproduced once the two-condition code is deleted. Re-running this script
# writes to `current_two_condition.tsv` instead, and refuses to overwrite the baseline.

using CarbonateCalculator

const BASELINE = "baseline_two_condition.tsv"
const OUT = joinpath(@__DIR__, isempty(ARGS) ? "current_two_condition.tsv" : ARGS[1])

if basename(OUT) == BASELINE && isfile(OUT)
    error("Refusing to overwrite $BASELINE - it is the pre-refactor reference and cannot " *
          "be regenerated. Write somewhere else, or delete it deliberately first.")
end

# name => (function, kwargs)
cases = Pair{String,Any}[]

push!(cases, "csys_T_20to2"        => (carbon_system, (TA=2300.0, DIC=2000.0, T_in=20.0, T_out=2.0)))
push!(cases, "csys_P_0to400"       => (carbon_system, (TA=2300.0, DIC=2000.0, T_in=20.0, P_in=0.0, P_out=400.0)))
push!(cases, "csys_TP_both"        => (carbon_system, (TA=2300.0, DIC=2000.0, T_in=20.0, T_out=2.0, P_in=0.0, P_out=400.0)))
push!(cases, "csys_S_28to38"       => (carbon_system, (TA=2300.0, DIC=2000.0, S_in=28.2, S_out=38.1)))
push!(cases, "csys_pHtot_TA"       => (carbon_system, (pHtot=8.1, TA=2300.0, T_in=25.0, T_out=5.0)))
push!(cases, "csys_pCO2_TA"        => (carbon_system, (pCO₂=400.0, TA=2300.0, T_in=25.0, T_out=5.0)))
push!(cases, "csys_nutrients"      => (carbon_system, (TA=2300.0, DIC=2000.0, PT=1.0, SiT=15.0, T_in=25.0, T_out=5.0)))
# bug case: H2ST/NH4T are not forwarded to the out-condition call
push!(cases, "csys_H2S_NH4"        => (carbon_system, (TA=2300.0, DIC=2000.0, H2ST=5.0, NH4T=3.0, T_in=25.0, T_out=5.0)))
# bug case: KSO4_method is hardcoded to "default" in the recursion
push!(cases, "csys_lueker_khoo"    => (carbon_system, (TA=2300.0, DIC=2000.0, T_in=20.0, T_out=2.0, K_method="Lueker 2000", KSO4_method="Khoo")))
push!(cases, "csys_lueker_default" => (carbon_system, (TA=2300.0, DIC=2000.0, T_in=20.0, T_out=2.0, K_method="Lueker 2000")))
push!(cases, "csys_unit_mmol"      => (carbon_system, (TA=2.3, DIC=2.0, unit="mmol", T_in=20.0, T_out=2.0)))
push!(cases, "csys_paleo_CaMg"     => (carbon_system, (TA=2300.0, DIC=2000.0, Ca=0.02, Mg=0.03, T_in=20.0, T_out=2.0)))

# bug case: δBT is not forwarded, so δBT_out comes back as the 39.61 default
# δBT is the modern seawater value, Isotopes.get_δBT() == 39.61.
const δBT_SW = 39.61

push!(cases, "wsys_dBT_T_20to2"    => (whole_system, (TA=2300.0, DIC=2000.0, δBT=δBT_SW, T_in=20.0, T_out=2.0)))
push!(cases, "wsys_pHtot_dBT"      => (whole_system, (pHtot=8.1, TA=2300.0, δBT=δBT_SW, T_in=20.0, T_out=2.0)))
push!(cases, "wsys_TP_both"        => (whole_system, (TA=2300.0, DIC=2000.0, δBT=δBT_SW, T_in=20.0, T_out=2.0, P_in=0.0, P_out=400.0)))
push!(cases, "wsys_S_28to38"       => (whole_system, (pHtot=8.1, TA=2300.0, δBT=δBT_SW, S_in=28.2, S_out=38.1)))
push!(cases, "wsys_alphaB"         => (whole_system, (TA=2300.0, DIC=2000.0, δBT=δBT_SW, alphaB=1.0272, T_in=20.0, T_out=2.0)))
# Deliberately NOT the modern value. δBT must be conserved across a condition change, and a
# case using the default cannot tell "carried correctly" apart from "silently re-defaulted"
# - which is exactly how the old recursion's dropped-δBT bug stayed hidden. Palaeo seawater
# genuinely differs, so this is a real configuration, not only a probe.
push!(cases, "wsys_dBT_nondefault" => (whole_system, (TA=2300.0, DIC=2000.0, δBT=38.0, T_in=20.0, T_out=2.0)))

open(OUT, "w") do io
    println(io, "case\tkey\tvalue")
    for (name, (f, kw)) in cases
        local res
        try
            res = f(; kw...)
        catch e
            println(io, "$name\t__ERROR__\t$(sprint(showerror, e))")
            @warn "case $name failed" exception = e
            continue
        end
        v = res.val
        for k in sort(collect(keys(v)); by = string)
            x = v[k]
            if x isa Number && isfinite(x)
                println(io, "$name\t$k\t", repr(Float64(x)))
            elseif x isa NamedTuple   # the Ks bundle
                for kk in sort(collect(keys(x)); by = string)
                    xx = x[kk]
                    xx isa Number && isfinite(xx) && println(io, "$name\t$k.$kk\t", repr(Float64(xx)))
                end
            end
        end
    end
end

n = count(==('\n'), read(OUT, String)) - 1
println("wrote $n rows to $OUT")
