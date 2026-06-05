using CSV, DataFrames
include("../../src/Calculator.jl")
using Calculator

# 1. K1 & K2 Constants
k1k2_map = Dict(
    1 => "Roy 1993", 2 => "GP 1989", 3 => "Hansson 1973", 4 => "HM 1973", 5 => "DM 1987",
    6 => "Mehrbach 1973 A", 7 => "Mehrbach 1973 B", 8 => "Millero 1979", 9 => "CW 2003", 10 => "Lueker 2000",
    11 => "MPM 2002", 12 => "Millero 2002", 13 => "Millero 2006", 14 => "Millero 2010", 
    15 => "Waters 2014", 16 => "Sulpis 2020", 17 => "SB 2020"
)

# 2. KSO4 Constants
kso4_map = Dict(1 => "Dickson", 2 => "Khoo", 3 => "Waters")

# 3. Total Boron Constants
boron_map = Dict(1 => "Uppstrom", 2 => "Lee")

# 4. Input Parameter Types
par_map = Dict(
    1 => :TA, 2 => :DIC, 3 => :pH_temp, 4 => :pCO₂, 
    5 => :fCO₂, 6 => :CO₃, 7 => :HCO₃, 8 => :CO₂
)

println("Starting test run across all rows...")

# Initialize the DataFrame to store failures with separated method columns
error_log = DataFrame(
    Row = Int[], 
    K_Method = String[], 
    KSO4_Method = String[], 
    B_Method = String[], 
    Parameter = String[], 
    Expected = Float64[], 
    Result = Float64[], 
    Difference = Float64[]
)

constants_log = DataFrame(
    Row = Int[], K_Method = String[],
    K0_exp = Float64[], K0_res = Float64[], K0_diff = Float64[],
    K1_exp = Float64[], K1_res = Float64[], K1_diff = Float64[],
    K2_exp = Float64[], K2_res = Float64[], K2_diff = Float64[],
    KW_exp = Float64[], KW_res = Float64[], KW_diff = Float64[],
    KB_exp = Float64[], KB_res = Float64[], KB_diff = Float64[],
    KF_exp = Float64[], KF_res = Float64[], KF_diff = Float64[],
    KS_exp = Float64[], KS_res = Float64[], KS_diff = Float64[],
    KP1_exp = Float64[], KP1_res = Float64[], KP1_diff = Float64[],
    KP2_exp = Float64[], KP2_res = Float64[], KP2_diff = Float64[],
    KP3_exp = Float64[], KP3_res = Float64[], KP3_diff = Float64[],
    KSi_exp = Float64[], KSi_res = Float64[], KSi_diff = Float64[],
    KH2S_exp = Float64[], KH2S_res = Float64[], KH2S_diff = Float64[],
    KNH3_exp = Float64[], KNH3_res = Float64[], KNH3_diff = Float64[]
)

k_rel_threshold = 1e-4 

k_anomaly_log = DataFrame(
    Row = Int[], 
    K1K2_Method = String[], 
    KSO4_Method = String[], 
    B_Method = String[], 
    K_Name = String[], 
    MATLAB_K = Float64[], 
    Julia_K = Float64[],
    Abs_Difference = Float64[],
    Rel_Difference = Float64[], # Now tracking relative difference
    T_in = Float64[], 
    S_in = Float64[], 
    P_in = Float64[]
)


# Dictionary to track Pass/Fail rates. Format: Combo_String => [Total_Tested, Total_Passed]
method_stats = Dict{String, Vector{Int}}()

df = CSV.read("results/compare_MATLABv3_2_0.csv", DataFrame)

for (i, row) in enumerate(eachrow(df))
    k_method = get(k1k2_map, row.K1K2CONSTANTS, missing)
    kso4_method = get(kso4_map, row.KSO4CONSTANT, missing)
    b_method = get(boron_map, row.BORON, missing)

    if ismissing(k_method) || ismissing(kso4_method) || ismissing(b_method)
        continue
    end

    # Keep a combo string specifically for the summary dictionary keys
    combo_name = "$k_method | $kso4_method | $b_method"
    
    # Initialize stats for this combination if it hasn't been seen yet
    if !haskey(method_stats, combo_name)
        method_stats[combo_name] = [0, 0]
    end
    
    # Increment total rows checked for this method combo
    method_stats[combo_name][1] += 1 

    inputs = Dict{Symbol, Any}()
    sym1 = par_map[row.PAR1TYPE]
    sym2 = par_map[row.PAR2TYPE]
    
    inputs[sym1] = row.PAR1
    inputs[sym2] = row.PAR2

    ph_scale = row.pHSCALEIN
    if haskey(inputs, :pH_temp)
        val = pop!(inputs, :pH_temp)
        if ph_scale == 1; inputs[:pHtot] = val; end
        if ph_scale == 2; inputs[:pHsws] = val; end
        if ph_scale == 3; inputs[:pHfree] = val; end
        if ph_scale == 4; inputs[:pHNBS] = val; end
    end

    env = (
        T_in = row.TEMPIN,
        S_in = row.SAL,
        P_in = row.PRESIN / 10.0, 
        PT = row.PO4,
        SiT = row.SI,
        H2ST = row.H2S,
        NH4T = row.NH4,
        unit = "umol",
        K_method = k_method, 
        KSO4_method = kso4_method,
        BT_method = b_method
    )

    try
        res = whole_system(; inputs..., env...)
        
        # Create a list of the K values we want to check
        k_checks = [
            ("K0", row.K0input, res.Ks.K0),
            ("K1", row.K1input, res.Ks.K1),
            ("K2", row.K2input, res.Ks.K2),
            ("KW", row.KWinput, res.Ks.KW),
            ("KB", row.KBinput, res.Ks.KB),
            ("KF", row.KFinput, res.Ks.KF),
            ("KS", row.KSinput, res.Ks.KS),
            ("KP1", row.KP1input, res.Ks.KP1),
            ("KP2", row.KP2input, res.Ks.KP2),
            ("KP3", row.KP3input, res.Ks.KP3),
            ("KSi", row.KSiinput, res.Ks.KSi),
            ("KH2S", row.KH2Sinput, res.Ks.KH2S),
            ("KNH3", row.KNH3input, res.Ks.KNH3)
        ]


        # Loop through them and log any that exceed the threshold
        for (k_name, exp_val, res_val) in k_checks
            # Calculate both absolute and relative differences
            abs_diff = res_val - exp_val
            rel_diff = exp_val != 0.0 ? abs(abs_diff / exp_val) : NaN
            
            # Trigger logging if the RELATIVE difference exceeds the threshold
            if !isnan(rel_diff) && rel_diff > k_rel_threshold
                push!(k_anomaly_log, (
                    i, k_method, kso4_method, b_method, k_name, 
                    exp_val, res_val, abs_diff, rel_diff, 
                    env.T_in, env.S_in, env.P_in
                ))
            end
        end


        row_passed = true
        
        # 1. Check DIC
        if !haskey(inputs, :DIC)
            expected = row.TCO2
            result = res.DIC
            if !isapprox(result, expected, atol=0.5)
                row_passed = false
                push!(error_log, (i, k_method, kso4_method, b_method, "DIC", expected, result, result - expected))
            end
        end
        
        # 2. Check pHtot
        if !haskey(inputs, :pHtot) && !haskey(inputs, :pHsws) && !haskey(inputs, :pHfree) && !haskey(inputs, :pHNBS)
            expected = row.pHinTOTAL
            result = res.pHtot
            if !isapprox(result, expected, atol=0.005)
                row_passed = false
                push!(error_log, (i, k_method, kso4_method, b_method, "pHtot", expected, result, result - expected))
            end
        end
        
        # 3. Check pCO2
        if !haskey(inputs, :pCO₂)
            expected = row.pCO2in
            result = res.pCO₂
            if !isapprox(result, expected, atol=0.5)
                row_passed = false
                push!(error_log, (i, k_method, kso4_method, b_method, "pCO2", expected, result, result - expected))
            end
        end
        
        # 4. Check CO3
        if !haskey(inputs, :CO₃)
            expected = row.CO3in
            result = res.CO₃
            if !isapprox(result, expected, atol=0.5)
                row_passed = false
                push!(error_log, (i, k_method, kso4_method, b_method, "CO3", expected, result, result - expected))
            end
        end
        
        if row_passed
            method_stats[combo_name][2] += 1 
        else
            # Log constants only for failed rows
            push!(constants_log, (
                i, k_method,
                row.K0input, res.Ks.K0, res.Ks.K0 - row.K0input,
                row.K1input, res.Ks.K1, res.Ks.K1 - row.K1input,
                row.K2input, res.Ks.K2, res.Ks.K2 - row.K2input,
                row.KWinput, res.Ks.KW, res.Ks.KW - row.KWinput,
                row.KBinput, res.Ks.KB, res.Ks.KB - row.KBinput,
                row.KFinput, res.Ks.KF, res.Ks.KF - row.KFinput,
                row.KSinput, res.Ks.KS, res.Ks.KS - row.KSinput,
                row.KP1input, res.Ks.KP1, res.Ks.KP1 - row.KP1input,
                row.KP2input, res.Ks.KP2, res.Ks.KP2 - row.KP2input,
                row.KP3input, res.Ks.KP3, res.Ks.KP3 - row.KP3input,
                row.KSiinput, res.Ks.KSi, res.Ks.KSi - row.KSiinput,
                row.KH2Sinput, res.Ks.KH2S, res.Ks.KH2S - row.KH2Sinput,
                row.KNH3input, res.Ks.KNH3, res.Ks.KNH3 - row.KNH3input
            ))
        end

    catch e
        # If the solver crashes, log it as a failure so it doesn't break the whole loop
        push!(error_log, (i, k_method, kso4_method, b_method, "CRASH", NaN, NaN, NaN))
    end
end

# --- FINALIZE AND EXPORT ---
println("Test run complete. Saving error log...")
CSV.write("error_log.csv", error_log)
CSV.write("error_log_constants.csv", constants_log)
CSV.write("error_log_k_anomalies.csv", k_anomaly_log)

println("\n--- RESULTS SUMMARY ---")
any_100_percent = false

for (combo, stats) in method_stats
    total_tested = stats[1]
    total_passed = stats[2]
    
    if total_tested > 0 && total_tested == total_passed
        println("[$combo] 100% PASS RATE")
        global any_100_percent = true
    end
end

if !any_100_percent
    println("No methods completely passed")
end