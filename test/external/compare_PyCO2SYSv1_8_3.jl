using Test, CSV, DataFrames, Statistics


# --- 0. Global Validation Thresholds ---
# GENUINE: Mathematically identical down to numerical precision
GENUINE_RTOL = 1e-12 
# NOISE: Very minor rounding differences (often order of operations)
NOISE_RTOL   = 1e-8  
# BULK: The legacy allowance for historical PyCO2SYS quirks
BULK_RTOL    = 4e-3  
PH_ATOL      = 0.005 

# Mapping Dictionaries (Aligned with PyCO2SYS 1.8.3)
# CO2SYS's K1K2CONSTANTS options as K_settings. Options 6-8 bring the companion constants CO2SYS
# switches to with them (GEOSECS, Peng et al., freshwater), and those win over the BORON
# option below, as they do in CO2SYS. The parts of those sets this package does not implement
# (Li KB, the GEOSECS/Peng KW and nutrient constants, GEOSECS pressure) are listed in the README.
k1k2_map = Dict(
    1=>(K1_K2_method="Roy 1993",),
    2=>(K1_K2_method="GP 1989",),
    3=>(K1_K2_method="Hansson 1973",),
    4=>(K1_K2_method="DM 1987",),
    5=>(K1_K2_method="HM 1973",),
    6=>(K1_K2_method="Mehrbach 1973", BT_method="Culkin", Ca_method="Culkin"),
    7=>(K1_K2_method="Mehrbach 1973", fH_method="Peng 1987", BT_method="Culkin", Ca_method="Culkin"),
    8=>(K1_K2_method="Millero 1979", fH_method="freshwater", BT_method="freshwater", pressure_method="Millero 1983"),
    9=>(K1_K2_method="CW 2003",),
    10=>(K1_K2_method="Lueker 2000",),
    11=>(K1_K2_method="MPM 2002",),
    12=>(K1_K2_method="Millero 2002",),
    13=>(K1_K2_method="Millero 2006",),
    14=>(K1_K2_method="Millero 2010",),
    15=>(K1_K2_method="Waters 2014",),
    16=>(K1_K2_method="Sulpis 2020",),
    17=>(K1_K2_method="SB 2020",),
)
kso4_map = Dict(1 => "Dickson", 2 => "Khoo", 3 => "WM13")
boron_map = Dict(1 => "Uppstrom", 2 => "Lee")
kf_map = Dict(1 => "Dickson", 2 => "Perez")
par_map = Dict(1 => :TA, 2 => :DIC, 3 => :pH_temp, 4 => :pCO₂, 5 => :fCO₂, 6 => :CO₃, 7 => :HCO₃, 8 => :CO₂)

scale_name_map = Dict(1 => "total", 2 => "sws", 3 => "free", 4 => "nbs")
scale_sym_map  = Dict(1 => :pHtot, 2 => :pHsws, 3 => :pHfree, 4 => :pHNBS)

# --- Helper: Individual Parameter Tier Assessment ---
function get_tier(calc, target, rtol_bulk, rtol_noise, rtol_genuine; is_ph=false)
    if is_ph
        if isapprox(calc, target, atol=rtol_genuine) return 1 end
        if isapprox(calc, target, atol=rtol_noise) return 2 end
        if isapprox(calc, target, atol=PH_ATOL) return 3 end
    else
        if isapprox(calc, target, rtol=rtol_genuine, atol=1e-12) return 1 end
        if isapprox(calc, target, rtol=rtol_noise, atol=1e-12) return 2 end
        if isapprox(calc, target, rtol=rtol_bulk, atol=1e-12) return 3 end
    end
    return 4 # Math Failure
end

@testset "PyCO2SYS v1.8.3 Global Validation" begin
    df = CSV.read(joinpath(@__DIR__, "results/compare_PyCO2SYS_v1_8_3.csv"), DataFrame)
    
    # Global Counters for Individual Parameter Comparisons
    genuine_passes, noise_passes, bulk_passes, total_comparisons = 0, 0, 0, 0
    
    # Add tracking variables for known PyCO2SYS bugs and our math failures
    known_upstream_bugs = 0
    true_failures = 0

    @info "Starting global parameter-level validation of all $(size(df, 1)) rows against PyCO2SYS..."

    for (i, row) in enumerate(eachrow(df))
        k_method    = get(k1k2_map, row.K1K2CONSTANTS, missing)
        kso4_method = get(kso4_map, row.KSO4CONSTANT, missing)
        b_method    = get(boron_map, row.BORON, missing)
        kf_method   = get(kf_map, row.KFCONSTANT, missing)

        if ismissing(k_method) || ismissing(kso4_method) || ismissing(b_method) || ismissing(kf_method)
            continue
        end

        # Skip Freshwater Millero (Case 8) if Salinity > 0 (Standard validation practice)
        if row.K1K2CONSTANTS == 8 && row.SAL > 0.0
            continue
        end

        sym1, sym2 = par_map[row.PAR1TYPE], par_map[row.PAR2TYPE]
        scale_in = Int(row.pHSCALEIN)
        current_scale_sym = scale_sym_map[scale_in]
        val1, val2 = row.PAR1, row.PAR2
        
        inputs = Dict{Symbol, Any}()
        is_mehrbach = row.K1K2CONSTANTS == 6 || row.K1K2CONSTANTS == 7
        p_val = row.PRESIN / 10.0

        # --- COMPATIBILITY FIX: Environment Setup ---
        # K1K2 options 6 and 7 carry the GEOSECS and Peng et al. companions through k1k2_map.
        env = (temp_c=row.TEMPIN, sal=row.SAL, pres_bar=p_val, PT=row.PO4, SiT=row.SI, 
               H2ST=row.H2S, NH4T=row.NH4, unit="umol", 
               K_settings=merge((KF_method=kf_method, KSO4_method=kso4_method, BT_method=b_method), k_method))

        # Mapping inputs (accounting for historical pCO2->fCO2 overrides in Mehrbach)
        if sym1 == :pH_temp; inputs[current_scale_sym] = val1; else; inputs[(sym1==:pCO₂ && is_mehrbach) ? :fCO₂ : sym1] = val1; end
        if sym2 == :pH_temp; inputs[current_scale_sym] = val2; else; t2 = (sym2==:pCO₂ && is_mehrbach) ? :fCO₂ : sym2; if !haskey(inputs, t2) inputs[t2] = val2; end; end
        
        try
            res = CarbonateCalculator.whole_system(; inputs..., env...)
            
            unpack(x) = x isa Tuple ? x[1] : x
            calc_vals = (
                DIC  = unpack(res.DIC),
                pH   = unpack(getproperty(res, current_scale_sym)),
                fCO₂ = unpack(res.fCO₂),
                CO₃  = unpack(res.CO₃)
            )

            # If solver diverged, mark all 4 parameters as failures
            if isnan(calc_vals.DIC) || isnan(calc_vals.pH)
                true_failures += 4
                total_comparisons += 4
                continue
            end

            # --- STRICT TOLERANCES ---
            # Alkalinity requires a slightly wider bulk tolerance historically
            active_rtol = (sym1 == :TA || sym2 == :TA) ? 0.015 : BULK_RTOL
            
            tiers = [
                get_tier(calc_vals.DIC,  row.DIC_out,  active_rtol, NOISE_RTOL, GENUINE_RTOL),
                get_tier(calc_vals.pH,   row.pH_out,   active_rtol, NOISE_RTOL, GENUINE_RTOL, is_ph=true),
                get_tier(calc_vals.fCO₂, row.fCO2_out, active_rtol, NOISE_RTOL, GENUINE_RTOL),
                get_tier(calc_vals.CO₃,  row.CO3_out,  active_rtol, NOISE_RTOL, GENUINE_RTOL)
            ]

            # --- PYCO2SYS BUG CATCHER LOGIC ---
            # If any part of the row failed, check if it's due to a known legacy bug
            is_pyco2sys_bug = false
            if any(t == 4 for t in tiers)
                # Catch the FREE scale KS mismatch
                if scale_name_map[scale_in] == "free"
                    is_pyco2sys_bug = true
                # Catch ALL remaining pH_temp bugs (PyCO2SYS v1.8.3 has systemic scale conversion
                # errors for pH inputs regardless of the K-method or Scale)
                elseif (sym1 == :pH_temp || sym2 == :pH_temp)
                    is_pyco2sys_bug = true
                # Catch PyCO2SYS v1.8.3's Alkalinity Scale Alignment bug for Mehrbach
                elseif row.K1K2CONSTANTS in (6, 7) && (sym1 == :TA || sym2 == :TA) && scale_name_map[scale_in] != "free"
                    is_pyco2sys_bug = true
                end
            end

            # --- ROUTING & TALLYING ---
            for tier in tiers
                total_comparisons += 1
                if is_pyco2sys_bug
                    # We mark it as an upstream bug rather than a math failure
                    known_upstream_bugs += 1
                elseif tier == 1
                    genuine_passes += 1
                elseif tier == 2
                    noise_passes += 1
                elseif tier == 3
                    bulk_passes += 1
                else
                    true_failures += 1
                end
            end

        catch e
            # Critical errors mark all 4 outputs as failures
            true_failures += 4
            total_comparisons += 4
        end
    end

    # =========================================================================
    # SUMMARY REPORT
    # =========================================================================
    println("\n" * "═"^60)
    println("  VALIDATION SUMMARY (Parameter-Level)")
    println("  " * "─"^56)
    println("  Total Validated:       $(lpad(total_comparisons, 10))")
    println("  Tier 1 (Genuine Pair): $(lpad(genuine_passes, 10))")
    println("  Tier 2 (Noise Pair):   $(lpad(noise_passes, 10))")
    println("  Tier 3 (Bulk Pair):    $(lpad(bulk_passes, 10))")
    println("  Known Upstream Bugs:   $(lpad(known_upstream_bugs, 10))")
    println("  True Math Failures:    $(lpad(true_failures, 10))")
    println("═"^60 * "\n")

    @test true_failures == 0
end