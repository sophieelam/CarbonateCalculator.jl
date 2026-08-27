using Printf

# `CarbonateResult` itself is defined in result.jl; this file only defines how it prints.

"""
Print a `CarbonateResult` as a summary block: the conditions, the parameters supplied, and
the quantities the calculation produced.

Inputs are left out of the body, so what is shown is what was computed rather than what was
handed in. The full state is always on `.val`.
"""
function Base.show(io::IO, ::MIME"text/plain", r::CarbonateResult)
    println(io, "══════════════════════════════════════════════════")
    println(io, "  CARBONATE SYSTEM RESULTS")
    println(io, "──────────────────────────────────────────────────")

    input_str = join([string(k) for k in r.input_keys], ", ")
    println(io, "  [ Inputs provided: ", input_str, " ]")
    @printf(io, "  [ Conditions: %.4g °C, S %.4g, %.4g bar ]\n",
            r.val.temp_c, r.val.sal, r.val.pres_bar)
    println(io, "──────────────────────────────────────────────────")

    _print_dynamic_vars(io, r)

    println(io, "──────────────────────────────────────────────────")
    println(io, "  [ Full results: .val | Uncertainties: .err ]")
    println(io, "══════════════════════════════════════════════════")
end

"""
Print the computed quantities of `r`, one line each, skipping anything the caller supplied.

Each row is looked up under several possible names so that the display survives a result
whose keys differ by scope. Precision follows the reporting unit: mol/kg needs ten decimal
places to show anything, µmol/kg needs four.
"""
function _print_dynamic_vars(io, r)
    # Define categories and potential keys found in the NamedTuple
    mapping = [
        "TA"             => [:TA, :Alk, :TAlk, :TotalAlkalinity],
        "DIC"            => [:DIC, :TC, :TCO2, :TotalCarbon],
        "pHtot"          => [:pHtot, :pH, :pH_tot],
        "pCO₂"           => [:pCO2, :pCO₂, :fCO2, :fCO₂],
        "CO₃"            => [:CO3, :CO₃],
        "HCO₃"           => [:HCO3, :HCO₃],
        "BOH₃"           => [:BOH3, :BOH₃, :BoricAcid],
        "BOH₄"           => [:BOH4, :BOH₄, :Borate],
        "δBOH₃"          => [:δBOH₃, :ABOH₃],
        "δBOH₄"          => [:δBOH₄, :ABOH₄],
        "ΩA"             => [:OmegaA, :ΩA, :OmegaAragonite],
        "ΩC"             => [:OmegaC, :ΩC, :OmegaCalcite],
    ]
    
    # --- Unit & Precision Handling ---
    unit_label = (r.unit == "umol" ? "μmol" : r.unit) * "/kg"

    # Every unit prints at the resolution µmol/kg gets at four decimal places, so the digits
    # shown follow the size of the numbers rather than the name they are reported under.
    decimals = max(round(Int, 4 + log10(1e6 / _unit_multiplier(r.unit))), 0)
    fmt_val = "%$(decimals + 4).$(decimals)f"
    fmt_err = "%.$(decimals)f"

    for (var_name, synonyms) in mapping
        # Skip if this variable category was an input
        if any(s -> s in r.input_keys, synonyms)
            continue
        end

        # Try to find a matching key
        found_key = nothing
        for s in synonyms
            if haskey(r.val, s)
                found_key = s
                break
            end
        end

        if found_key !== nothing
            val = r.val[found_key]
            
            # Skip boron variables if they are empty/zero (standard system calls)
            boron_names = ["BOH₃", "BOH₄", "δBOH₃", "δBOH₄"]
            if var_name in boron_names && (val === nothing || val == 0.0 || isnan(val))
                continue
            end

            err = (r.err !== nothing && haskey(r.err, found_key)) ? r.err[found_key] : 0.0
            
            # Construct the Pretty Label
            label = var_name == "pHtot"     ? "pH (Total)" :
                    var_name == "pCO₂"      ? "pCO₂ (μatm)" :
                    var_name == "CO₃"       ? "CO₃ ($unit_label)" :
                    var_name == "HCO₃"      ? "HCO₃ ($unit_label)" :
                    var_name == "TA"        ? "TA ($unit_label)" :
                    var_name == "DIC"       ? "DIC ($unit_label)" :
                    var_name == "BOH₃"      ? "B(OH)₃ ($unit_label)" :
                    var_name == "BOH₄"      ? "B(OH)₄⁻ ($unit_label)" :
                    var_name == "δBOH₃"     ? "δ¹¹B_B(OH)₃ (‰)" :
                    var_name == "δBOH₄"     ? "δ¹¹B_B(OH)₄⁻ (‰)" :
                    var_name == "ΩA"        ? "Ω Aragonite" :
                    var_name == "ΩC"        ? "Ω Calcite" : "Unknown"
            
            label_padded = rpad(label, 22)
            
            if err > 0.0
                # Specific formatting for isotopes
                if var_name in ("δBOH₃", "δBOH₄")
                    @printf(io, "    %s : %8.2f ± %.3f\n", label_padded, val, err)
                else
                    # Dynamic precision based on unit scale
                    val_str = Printf.format(Printf.Format("    %s : $fmt_val ± $fmt_err\n"), label_padded, val, err)
                    print(io, val_str)
                end
            else
                val_str = Printf.format(Printf.Format("    %s : $fmt_val\n"), label_padded, val)
                print(io, val_str)
            end
        end
    end
end