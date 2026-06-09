using Printf

# 1. Define the custom struct
struct CarbonateResult
    val::NamedTuple
    err::Union{NamedTuple, Nothing}
    input_keys::Vector{Symbol} 
end

function Base.getproperty(res::CarbonateResult, s::Symbol)
    # 1. Standard fields
    if s in (:val, :err, :input_keys)
        return getfield(res, s)
    end
    
    v = getfield(res, :val)
    
    # 2. Priority: If it's a two-state system, return the "out" value as the default
    s_out = Symbol(string(s), "_out")
    if haskey(v, s_out)
        return getproperty(v, s_out)
    end
    
    # 3. Fallback: Return the base value (or the "in" value if that's all there is)
    return getproperty(v, s)
end

# Make tab-completion work for both the struct fields AND the math outputs
function Base.propertynames(res::CarbonateResult, private::Bool=false)
    return (fieldnames(CarbonateResult)..., keys(getfield(res, :val))...)
end

# Allow iteration so tools like ForwardDiff can treat the result like a tuple
function Base.iterate(res::CarbonateResult, state...)
    return iterate(getfield(res, :val), state...)
end

function Base.keys(res::CarbonateResult)
    return keys(getfield(res, :val))
end

function Base.length(res::CarbonateResult)
    return length(getfield(res, :val))
end

# 2. Overload the Base.show method
function Base.show(io::IO, ::MIME"text/plain", r::CarbonateResult)
    is_two_state = haskey(r.val, :pHtot_in)
    
    println(io, "══════════════════════════════════════════════════")
    println(io, "  CARBONATE SYSTEM RESULTS (Dynamic View)")
    println(io, "──────────────────────────────────────────────────")
    
    input_str = join([string(k) for k in r.input_keys], ", ")
    println(io, "  [ Inputs provided: ", input_str, " ]")
    println(io, "──────────────────────────────────────────────────")

    if is_two_state
        println(io, "  [ Surface / Input State ]")
        _print_dynamic_vars(io, r, "_in")
        println(io, "──────────────────────────────────────────────────")
        println(io, "  [ Output State ]")
        _print_dynamic_vars(io, r, "_out")
    else
        _print_dynamic_vars(io, r, "")
    end

    println(io, "──────────────────────────────────────────────────")
    println(io, "  [ Full results: .val | Uncertainties: .err ]")
    println(io, "══════════════════════════════════════════════════")
end

function _print_dynamic_vars(io, r, suffix)
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
        "deltaBOH₄"      => [:deltaBOH4, :deltaBOH₄, :d11B_BOH4],
        "ΩA"             => [:OmegaA, :ΩA, :OmegaAragonite],
        "ΩC"             => [:OmegaC, :ΩC, :OmegaCalcite],
        "revelle_factor" => [:revelle_factor, :Revelle, :RF]
    ]
    
    # --- Unit & Precision Handling ---
    # --- Unit & Precision Handling ---
    # We check if :unit exists, and convert it to string to handle Symbols or Strings
    raw_unit = haskey(r.val, :unit) ? string(r.val[:unit]) : "umol"
    
    if raw_unit == "mol" || raw_unit == "mol/kg"
        unit_label = "mol/kg"
        fmt_val = "%14.10f"  # Extra width and 10 decimals for molar
        fmt_err = "%.10f"
    else
        unit_label = "μmol/kg"
        fmt_val = "%8.4f"    # 4 decimals for micromolar
        fmt_err = "%.4f"
    end

    for (var_name, synonyms) in mapping
        # Skip if this variable category was an input
        if any(s -> s in r.input_keys, synonyms)
            continue
        end

        # Try to find a matching key
        found_key = nothing
        for s in synonyms
            s_suffix = Symbol(string(s) * suffix)
            
            # 1. Look for the suffixed version (e.g., DIC_out)
            # 1. Look for the suffixed version (e.g., DIC_out)
            if haskey(r.val, s_suffix)
                found_key = s_suffix
                break
            # 2. FALLBACK:
            elseif suffix == "_out"
                if haskey(r.val, s)        # Look for the base symbol first
                    found_key = s
                    break
                else 
                    s_in = Symbol(string(s) * "_in")
                    if haskey(r.val, s_in) # Then fall back to _in
                        found_key = s_in
                        break
                    end
                end
            end
        end

        if found_key !== nothing
            val = r.val[found_key]
            
            # Skip boron variables if they are empty/zero (standard system calls)
            boron_names = ["BOH₃", "BOH₄", "deltaBOH₄"]
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
                    var_name == "deltaBOH₄" ? "δ¹¹B_borate (‰)" :
                    var_name == "ΩA"        ? "Ω Aragonite" :
                    var_name == "ΩC"        ? "Ω Calcite" : 
                    var_name == "revelle_factor" ? "Revelle Fact." : var_name
            
            label_padded = rpad(label, 22)
            
            if err > 0.0
                # Specific formatting for isotopes
                if var_name == "deltaBOH₄"
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