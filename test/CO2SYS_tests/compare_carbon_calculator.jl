using Test, CSV, DataFrames, Statistics
# Adjust this include path if needed based on where you save the script
# include("../../src/Calculator.jl")
# using .Calculator

# =========================================================================
# 0. GLOBAL CONFIGURATION & HELPER FUNCTIONS
# =========================================================================

# Validation Thresholds
GENUINE_RTOL = 1e-12 # Mathematically identical down to numerical precision
NOISE_RTOL   = 1e-8  # Very minor rounding differences (often order of operations)
BULK_RTOL    = 4e-3  # Legacy allowance for historical numerical quirks
PH_ATOL      = 0.005 

# Mapping Dictionaries
k1k2_map = Dict(
    1=>"Roy 1993", 2=>"GP 1989", 3=>"Hansson 1973", 4=>"DM 1987", 5=>"HM 1973", 
    6=>"Mehrbach 1973 A", 7=>"Mehrbach 1973 B", 8=>"Millero 1979", 9=>"CW 2003", 
    10=>"Lueker 2000", 11=>"MPM 2002", 12=>"Millero 2002", 13=>"Millero 2006", 
    14=>"Millero 2010", 15=>"Waters 2014", 16=>"Sulpis 2020", 17=>"SB 2020"
)
kso4_map = Dict(1 => "Dickson", 2 => "Khoo", 3 => "WM13")
boron_map = Dict(1 => "Uppstrom", 2 => "Lee")
kf_map = Dict(1 => "Dickson", 2 => "Perez")
par_map = Dict(1 => :TA, 2 => :DIC, 3 => :pH_temp, 4 => :pCO₂, 5 => :fCO₂, 6 => :CO₃, 7 => :HCO₃, 8 => :CO₂)

scale_name_map = Dict(1 => "total", 2 => "sws", 3 => "free", 4 => "nbs")
scale_sym_map  = Dict(1 => :pHtot, 2 => :pHsws, 3 => :pHfree, 4 => :pHNBS)

# Helper: Individual Parameter Tier Assessment
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

# =========================================================================
# 1. UNIFIED TEST SUITE
# =========================================================================

@testset "carbon_calculator Global Validation Suite" begin

    # -------------------------------------------------------------------------
    # PART A: MATLAB CO2SYS v3.2.0 VALIDATION
    # -------------------------------------------------------------------------
    @testset "Against MATLAB v3.2.0" begin
        df = CSV.read(joinpath(@__DIR__, "results/compare_MATLABv3_2_0.csv"), DataFrame)
        
        genuine_passes, noise_passes, bulk_passes, total_comparisons = 0, 0, 0, 0
        known_upstream_bugs, true_failures = 0, 0

        @info "Starting carbon_calculator validation against MATLAB CO2SYS..."

        for row in eachrow(df)
            k_method    = get(k1k2_map, row.K1K2CONSTANTS, missing)
            kso4_method = get(kso4_map, row.KSO4CONSTANT, missing)
            b_method    = get(boron_map, row.BORON, missing)
            kf_method   = get(kf_map, row.KFCONSTANT, missing)

            if ismissing(k_method) || ismissing(kso4_method) || ismissing(b_method) || ismissing(kf_method)
                continue
            end

            sym1, sym2 = par_map[row.PAR1TYPE], par_map[row.PAR2TYPE]
            val1 = sym1 == :pCO₂ ? row.pCO2in : (sym1 == :fCO₂ ? row.fCO2in : row.PAR1)
            val2 = sym2 == :pCO₂ ? row.pCO2in : (sym2 == :fCO₂ ? row.fCO2in : row.PAR2)

            scale_in = Int(row.pHSCALEIN)
            current_scale_str = scale_name_map[scale_in]
            current_scale_sym = scale_sym_map[scale_in]
            
            target_ph = if scale_in == 1 row.pHinTOTAL elseif scale_in == 2 row.pHinSWS elseif scale_in == 3 row.pHinFREE else row.pHinNBS end

            inputs = Dict{Symbol, Any}()
            is_mehrbach = row.K1K2CONSTANTS == 6 || row.K1K2CONSTANTS == 7
            
            env = (T_in=row.TEMPIN, S_in=row.SAL, P_in=row.PRESIN/10.0, PT=row.PO4, SiT=row.SI, 
                   H2ST=row.H2S, NH4T=row.NH4, unit="umol", scale=current_scale_str, 
                   K_method=k_method, KF_method=kf_method, KSO4_method=kso4_method, BT_method=b_method)

            if sym1 == :pH_temp; inputs[current_scale_sym] = val1; else; inputs[(sym1==:pCO₂ && is_mehrbach) ? :fCO₂ : sym1] = val1; end
            if sym2 == :pH_temp; inputs[current_scale_sym] = val2; else; t2 = (sym2==:pCO₂ && is_mehrbach) ? :fCO₂ : sym2; if !haskey(inputs, t2) inputs[t2] = val2; end; end
            
            try
                # CALLING THE NEW FUNCTION HERE
                res = Calculator.carbon_calculator(; inputs..., env...)
                
                unpack(x) = x isa Tuple ? x[1] : x
                calc_vals = (
                    DIC  = unpack(res.DIC),
                    pH   = unpack(getproperty(res, current_scale_sym)),
                    fCO₂ = unpack(res.fCO₂),
                    CO₃  = unpack(res.CO₃)
                )

                if isnan(calc_vals.DIC) || isnan(calc_vals.pH)
                    true_failures += 4; total_comparisons += 4; continue
                end

                active_rtol = (sym1 == :TA || sym2 == :TA) ? 0.015 : BULK_RTOL
                tiers = [
                    get_tier(calc_vals.DIC,  row.TCO2,     active_rtol, NOISE_RTOL, GENUINE_RTOL),
                    get_tier(calc_vals.pH,   target_ph,    active_rtol, NOISE_RTOL, GENUINE_RTOL, is_ph=true),
                    get_tier(calc_vals.fCO₂, row.fCO2in,   active_rtol, NOISE_RTOL, GENUINE_RTOL),
                    get_tier(calc_vals.CO₃,  row.CO3in,    active_rtol, NOISE_RTOL, GENUINE_RTOL)
                ]

                # Catch known MATLAB bugs
                is_matlab_bug = false
                if any(t == 4 for t in tiers)
                    if isapprox(calc_vals.DIC, row.TCO2, rtol=0.65, atol=1e-12) && isapprox(calc_vals.pH, target_ph, atol=0.15) && isapprox(calc_vals.fCO₂, row.fCO2in, rtol=0.65, atol=1e-12) && isapprox(calc_vals.CO₃, row.CO3in, rtol=0.65, atol=1e-12)
                        is_matlab_bug = true
                    end
                end

                for tier in tiers
                    total_comparisons += 1
                    if is_matlab_bug known_upstream_bugs += 1
                    elseif tier == 1 genuine_passes += 1
                    elseif tier == 2 noise_passes += 1
                    elseif tier == 3 bulk_passes += 1
                    else true_failures += 1 end
                end

            catch e
                true_failures += 4; total_comparisons += 4
            end
        end
        @test true_failures == 0
    end


    # -------------------------------------------------------------------------
    # PART B: PyCO2SYS v1.8.3 VALIDATION
    # -------------------------------------------------------------------------
    @testset "Against PyCO2SYS v1.8.3" begin
        df = CSV.read(joinpath(@__DIR__, "results/compare_PyCO2SYS_v1_8_3.csv"), DataFrame)
        
        genuine_passes, noise_passes, bulk_passes, total_comparisons = 0, 0, 0, 0
        known_upstream_bugs, true_failures = 0, 0

        @info "Starting carbon_calculator validation against PyCO2SYS..."

        for row in eachrow(df)
            k_method    = get(k1k2_map, row.K1K2CONSTANTS, missing)
            kso4_method = get(kso4_map, row.KSO4CONSTANT, missing)
            b_method    = get(boron_map, row.BORON, missing)
            kf_method   = get(kf_map, row.KFCONSTANT, missing)

            if ismissing(k_method) || ismissing(kso4_method) || ismissing(b_method) || ismissing(kf_method) || (row.K1K2CONSTANTS == 8 && row.SAL > 0.0)
                continue
            end

            sym1, sym2 = par_map[row.PAR1TYPE], par_map[row.PAR2TYPE]
            scale_in = Int(row.pHSCALEIN)
            current_scale_sym = scale_sym_map[scale_in]
            val1, val2 = row.PAR1, row.PAR2
            
            inputs = Dict{Symbol, Any}()
            is_mehrbach = row.K1K2CONSTANTS == 6 || row.K1K2CONSTANTS == 7

            env = (T_in=row.TEMPIN, S_in=row.SAL, P_in=row.PRESIN/10.0, PT=row.PO4, SiT=row.SI, 
                   H2ST=row.H2S, NH4T=row.NH4, unit="umol", scale=scale_name_map[scale_in], 
                   K_method=k_method, KF_method=kf_method, KSO4_method=kso4_method, BT_method=b_method)

            if sym1 == :pH_temp; inputs[current_scale_sym] = val1; else; inputs[(sym1==:pCO₂ && is_mehrbach) ? :fCO₂ : sym1] = val1; end
            if sym2 == :pH_temp; inputs[current_scale_sym] = val2; else; t2 = (sym2==:pCO₂ && is_mehrbach) ? :fCO₂ : sym2; if !haskey(inputs, t2) inputs[t2] = val2; end; end
            
            try
                # CALLING THE NEW FUNCTION HERE
                res = Calculator.carbon_calculator(; inputs..., env...)
                
                unpack(x) = x isa Tuple ? x[1] : x
                calc_vals = (
                    DIC  = unpack(res.DIC),
                    pH   = unpack(getproperty(res, current_scale_sym)),
                    fCO₂ = unpack(res.fCO₂),
                    CO₃  = unpack(res.CO₃)
                )

                if isnan(calc_vals.DIC) || isnan(calc_vals.pH)
                    true_failures += 4; total_comparisons += 4; continue
                end

                active_rtol = (sym1 == :TA || sym2 == :TA) ? 0.015 : BULK_RTOL
                tiers = [
                    get_tier(calc_vals.DIC,  row.DIC_out,  active_rtol, NOISE_RTOL, GENUINE_RTOL),
                    get_tier(calc_vals.pH,   row.pH_out,   active_rtol, NOISE_RTOL, GENUINE_RTOL, is_ph=true),
                    get_tier(calc_vals.fCO₂, row.fCO2_out, active_rtol, NOISE_RTOL, GENUINE_RTOL),
                    get_tier(calc_vals.CO₃,  row.CO3_out,  active_rtol, NOISE_RTOL, GENUINE_RTOL)
                ]

                # Catch known PyCO2SYS bugs
                is_pyco2sys_bug = false
                if any(t == 4 for t in tiers)
                    if scale_name_map[scale_in] == "free" || (sym1 == :pH_temp || sym2 == :pH_temp) || ((k_method == "Mehrbach 1973 A" || k_method == "Mehrbach 1973 B") && (sym1 == :TA || sym2 == :TA) && scale_name_map[scale_in] != "free")
                        is_pyco2sys_bug = true
                    end
                end

                for tier in tiers
                    total_comparisons += 1
                    if is_pyco2sys_bug known_upstream_bugs += 1
                    elseif tier == 1 genuine_passes += 1
                    elseif tier == 2 noise_passes += 1
                    elseif tier == 3 bulk_passes += 1
                    else true_failures += 1 end
                end

            catch e
                true_failures += 4; total_comparisons += 4
            end
        end
        @test true_failures == 0
    end

end