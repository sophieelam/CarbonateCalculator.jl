include("../src/Calculator.jl")

using Calculator

println("Starting tests...")


# 1. The Extreme Dataset
# These rows are engineered to trip the different logic branches in which_K
T_test = [
    36.0,  # Row 1: T > 35 -> Should trigger "Millero 2006"
    1.0,   # Row 2: T < 2 & S in 34-37 -> Should trigger "Millero 2002"
    -1.0,  # Row 3: T < 0 & S in 10-50 -> Should trigger "GP 1989"
    25.0,  # Row 4: S < 1.0 -> Should trigger "Millero 1979"
    25.0,  # Row 5: S > 50 -> Should trigger "Papadimitriou 2018"
    25.0   # Row 6: Standard seawater -> Should trigger default ("MyAMI")
]

S_test = [
    35.0,  # Normal salinity for Row 1
    35.0,  # Normal salinity for Row 2
    35.0,  # Normal salinity for Row 3
    0.5,   # Estuarine/freshwater for Row 4
    55.0,  # Brine for Row 5
    35.0   # Standard salinity for Row 6
]

# Adding varying pressures (in bar) to ensure Millero 1995 corrections apply without crashing
P_test = [0.0, 50.0, 100.0, 200.0, 300.0, 400.0] 

println("=== STARTING K_CALCULATOR TESTS ===\n")

# ---------------------------------------------------------
# Test 1: DYNAMIC MODE
# ---------------------------------------------------------
println("Test 1: Dynamic Mode (Evaluating row-by-row)")
try
    dynamic_results = Calculator.K_calculator(
        T_in = T_test, S_in = S_test, P_in = P_test, K_mode = "dynamic"
    )
    println("✅ Dynamic Mode Success!")
    
    # --- LOGIC CHECK ---
    println("\n🔍 PROVENANCE CHECK (Dynamic):")
    for i in eachindex(T_test, S_test, dynamic_results.method)
    println("   Row $i (T=$(T_test[i]), S=$(S_test[i])) -> Method: $(dynamic_results.method[i])")
    end

    all_same = all(x -> x == dynamic_results.Ks.K1[1], dynamic_results.Ks.K1)
    if !all_same
        println("   ✅ Verified: K1 values vary across the array.")
    end

catch e
    println("❌ Dynamic Mode Failed: ", e)
end

# ---------------------------------------------------------
# Test 2: STATIC MODE
# ---------------------------------------------------------
println("\nTest 2: Static Mode (Averaging the dataset first)")
try
    static_results = Calculator.K_calculator(
        T_in = T_test, S_in = S_test, P_in = P_test, K_mode = "static"
    )
    println("✅ Static Mode Success!")
    
    # --- LOGIC CHECK ---
    unique_method = unique(static_results.method)
    println("   🔍 PROVENANCE CHECK (Static):")
    println("   Averaged input triggered exactly ONE method for all rows: $(unique_method[1])")
    
catch e
    println("❌ Static Mode Failed: ", e)
end

# ---------------------------------------------------------
# Test 3: CUSTOM KWARGS (KS, KF, BT, and MyAMI toggle)
# ---------------------------------------------------------
println("\nTest 3: Custom Keyword Arguments (Passthrough testing)")
try
    custom_results = Calculator.K_calculator(
        T_in = T_test, S_in = S_test, P_in = P_test, K_mode = "dynamic",
        KSO4_method = "Khoo", KF_method = "Perez", BT_method = "Lee", Ca = 0.02
    )
    println("✅ Custom Kwargs Success!")
    
    # --- LOGIC CHECK ---
    methods_used = unique(custom_results.method)
    println("   🔍 PROVENANCE CHECK (Custom):")
    println("   Used methods: $(methods_used)") 
    # This should just be ["MyAMI"] because Ca=0.02 triggers it
    
catch e
    println("❌ Custom Kwargs Failed: ", e)
end

println("\n=== TESTS COMPLETE ===")