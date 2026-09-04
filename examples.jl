using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Revise
using CarbonateCalculator
using DataFrames
using Plots


## =============================================================================
# SECTION 0: TL;DR-- basic & advanced usage example walk-throughs
# ==============================================================================

# Basic usage example with keyword argumemnts & subsystem functions

res = carbon_system(TA=2300, DIC=2100)

df = DataFrame([res])

df.pHtot

# Advanced usage example with positional argumets and CarbonSystem

# first, build a solver for the specific subsystem(s) you want to solve for, and the parameters you intend to give it.
solver = CarbonateSystem(
    :isotopes,  # the subsystem(s) to solve for, in this case B isotopes
    varying=(:δBOH₄,),  # the parameters you will provide to the solver, in this case δ¹¹B of borate
    varying_errors=(:δBOH₄,),  # the parameters you will provide errors for, in this case δ¹¹B of borate
    sal=35.0, temp_c=20.0  # the conditions that are the same for every calculation, in this case salinity and temperature
    )

# this produces a callable solver, which accepts the `varying` and `varying_errors` parameters in the order that they were provided.
# This is a broadcastable function, which accepts arrays of values for the varying parameters, and returns an array of results.

result = solver.(  # broadcast the solver over the following values
    [15, 16, 17, 18],  # the values of δ¹¹B of borate
    [0.5, 0.5, 0.5, 0.5]  # the errors of δ¹¹B of borate, which is the same for every calculation in this case
    )

df = DataFrame(result)  # convert the results to a DataFrame for easier viewing and analysis

df[!, [:pHtot, :σ_pHtot]]  # view the calculated pH and its resulting uncertainty  

result[1]  # view a single result, accessed by index

## =============================================================================
# SECTION 1: Core Speciation Presets (Quick Keyword Usage)
# ==============================================================================
# Presets are pre-built solvers ideal for quick interactive calculations.
# Supply any two independent carbonate or isotopic constraints as keywords.

# 1a. Carbonate system only
res_carbon = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0, sal = 35.0)
println("pH Total (Carbon only): ", res_carbon.pHtot)

# 1b. Whole system (Carbon + Boron + Boron Isotopes)
res_whole = whole_system(TA = 2300.0, δBOH₄ = 16.0, temp_c = 20.0)
println("pH Total (Whole system from TA & δBOH₄): ", res_whole.pHtot)

# 1c. Boron speciation and isotopes only (no carbon machinery)
res_boron = boron_system(pHtot = 8.1, temp_c = 20.0)
println("B(OH)₄⁻ concentration: ", res_boron.BOH₄, " μmol/kg")

# 1d. Boron isotopes only
res_iso = boron_isotopes(pHtot = 8.1, temp_c = 20.0)
println("δ¹¹B of B(OH)₄⁻: ", res_iso.δBOH₄, " ‰")


## =============================================================================
# SECTION 2: Vectorized Solvers & DataFrames Integration
# ==============================================================================
# Keyword presets cannot broadcast over vectors. To solve multiple samples,
# construct a `CarbonateSystem` solver specifying `varying` input arguments.

# Create a compiled solver for TA, DIC, and Temperature
solver_vec = CarbonateSystem(
    :carbon;
    varying = (:TA, :DIC, :temp_c),
    K_method = "Lueker 2000",
    sal = 35.0
)

# Sample dataset
ta_vec   = [2280.0, 2300.0, 2320.0]
dic_vec  = [2000.0, 2050.0, 2100.0]
temp_vec = [25.0,   18.0,   10.0]

# Broadcast the solver over the input vectors
results_vec = solver_vec.(ta_vec, dic_vec, temp_vec)

# Convert results directly into a Tables.jl-compatible DataFrame
df_results = DataFrame(results_vec)
println("Vectorized DataFrame pH: ", df_results.pHtot)


## =============================================================================
# SECTION 3: Exact Uncertainty Propagation (Automatic Differentiation)
# ==============================================================================
# Track uncertainties without analytical formulas by naming uncertain parameters
# in `varying_errors`. ForwardDiff.jl evaluates exact partial derivatives.

solver_err = CarbonateSystem(
    :carbon;
    varying = (:TA, :DIC),
    varying_errors = (:TA, :DIC),
    temp_c = 25.0,
    sal = 35.0
)

# Arguments: (TA_val, DIC_val, σ_TA, σ_DIC)
res_err = solver_err(2300.0, 2000.0, 2.0, 2.0)

println("pH Total: ", res_err.pHtot)
println("pH Uncertainty (σ_pHtot): ", res_err.err.pHtot)


## =============================================================================
# SECTION 4: In-Situ & Collection Condition Adjustments
# ==============================================================================
# Samples are often measured at deck temperature/pressure (e.g., 25°C, 0 bar)
# but need to be reported at collection depth (in-situ).

measured_deck = carbon_system(TA = 2300.0, DIC = 2050.0, temp_c = 25.0, pres_bar = 0.0)

# Re-solve conservative totals at collection conditions (4°C, 250 bar)
in_situ = at_collection_conditions(measured_deck; temp_c = 4.0, pres_bar = 250.0)

println("Deck pH (25°C, 0 bar):   ", measured_deck.pHtot)
println("In-Situ pH (4°C, 250 bar): ", in_situ.pHtot)


# ==============================================================================
# SECTION 5: Chemical Gradients & Buffer Factors
# ==============================================================================
# Calculate thermodynamic buffer factors (e.g., Revelle factor) and partial derivatives.

res_buffer = carbon_system(TA = 2300.0, DIC = 2080.0, temp_c = 15.0)

# Revelle factor: ∂ln(pCO₂)/∂ln(DIC) at constant TA
rev_fac = revelle_factor(res_buffer)
println("Revelle Factor: ", rev_fac)

# Absolute gradient: ∂(pH)/∂(TA) at constant DIC
dpH_dTA_abs = calc_gradient(res_buffer, :pHtot, :TA, constant = :DIC)
println("∂(pH)/∂(TA) [constant DIC]: ", dpH_dTA_abs)

# Relative gradient: ∂(pH)/∂(ln TA) at constant DIC
dpH_dTA_rel = calc_relative_gradient(res_buffer, :pHtot, :TA, constant = :DIC)
println("∂(pH)/∂(ln TA) [constant DIC]: ", dpH_dTA_rel)

## =============================================================================
# SECTION 6: Calling different constant calculation methods
# ==============================================================================

# Explicitly define every parameterisation choice in a single call
res_custom_methods = carbon_system(
    TA = 2300.0, DIC = 2000.0, temp_c = 25.0, sal = 35.0,
    K_method    = "Lueker 2000",          # Carbonic acid constants
    BT_method   = "Lee",             # Total boron to salinity relationship
    KSO4_method = "Dickson",         # Sulphate dissociation constant
    KF_method   = "Perez", # Fluoride dissociation constant
    KNH3_method = "Clegg",           # Ammonia dissociation constant
    Ca_method   = "RT67",           # Calcium concentration relationship
    MyAMI_mode  = "calculate"             # Ca/Mg activity coefficient corrections
)

# Or pre-compute constants specifying all parameterisation methods:
ks_custom = CarbonateCalculator.calculate_constants(
    temp_c = 25.0, sal = 35.0, pres_bar = 0.0,
    K_method    = "Roy 1993",
    BT_method   = "Uppstrom",
    KSO4_method = "Khoo",
    KF_method   = "Dickson",
    KNH3_method = "Clegg",
    Ca_method   = "RT67",
    MyAMI_mode  = "approximate"
)

# Pass the pre-computed bundle directly using `Ks`
res_from_ks = carbon_system(TA = 2300.0, DIC = 2000.0, Ks = ks_custom)
println("pH using pre-computed custom Ks: ", res_from_ks.pHtot)


## =============================================================================
# SECTION 7: Paleoclimate Case Study
# ==============================================================================
# Scenario: A paleo-oceanographic hydrocast depth profile.
# We reconstruct seawater pCO₂ and pH down the water column from measured δ¹¹B
# of borate and TA, propagating analytical uncertainties, adjusting for in-situ
# collection temperature and pressure, and computing buffer factors.

# 1. Define Hydrocast Depth Profile
depths_m     = [10.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0]
pres_bars    = depths_m ./ 10.0                          # Hydrostatic pressure (bar)
deck_temp    = fill(25.0, length(depths_m))             # Deck measurement temperature (°C)
insitu_temp  = [22.0, 15.0, 8.0, 4.0, 2.5, 1.5]          # In-situ water temperature (°C)

# Measured variables on deck with analytical uncertainties (σ)
ta_meas      = [2290.0, 2300.0, 2310.0, 2325.0, 2340.0, 2350.0] # μmol/kg
d11b_meas    = [17.5, 17.1, 16.5, 15.8, 15.2, 15.0]           # ‰
σ_ta         = fill(2.0, length(depths_m))                   # ±2.0 μmol/kg
σ_d11b       = fill(0.2, length(depths_m))                   # ±0.2 ‰

# 2. Build Full Multi-Subsystem Solver with Error Propagation
mega_solver = CarbonateSystem(
    :carbon, :boron, :isotopes;                 # All three subsystems active
    varying = (:TA, :δBOH₄, :temp_c),             # Inputs provided positionally
    varying_errors = (:TA, :δBOH₄),             # Inputs with uncertainties
    sal = 35.0,
    K_method = "KGen"
)

# 3. Solve Deck Conditions with AD Uncertainty Propagation
deck_results = mega_solver.(ta_meas, d11b_meas, deck_temp, σ_ta, σ_d11b)

# 4. Project Results to In-Situ Collection Conditions (Depth, Pressure, Temp)
insitu_results = at_collection_conditions.(deck_results, insitu_temp, pres_bars)

# 5. Extract In-Situ Revelle Factors
revelle_insitu = revelle_factor.(insitu_results)

# 6. Build Master DataFrame
df_mega = DataFrame(insitu_results)
df_mega[!, :depth_m]       = depths_m
df_mega[!, :RevelleFactor] = revelle_insitu

# Select key oceanographic summary columns
summary_df = df_mega[!, [
    :depth_m, :pHtot, :σ_pHtot, :pCO₂, :σ_pCO₂, :ΩA, :σ_ΩA, :RevelleFactor
]]

println("\n=== MEGA BOSS HYDROCAST PROFILE RESULTS ===")
println(summary_df)

# 7. Plot Depth Profile with Analytical Uncertainty Ribbons
p1 = plot(
    df_mega.pHtot, df_mega.depth_m,
    xerror = df_mega.σ_pHtot,
    yflip = true,
    xlabel = "In-Situ pH (Total Scale)",
    ylabel = "Depth (m)",
    label = "pH ± σ",
    marker = :circle,
    linewidth = 2,
    color = :navy
)

p2 = plot(
    df_mega.pCO₂, df_mega.depth_m,
    xerror = df_mega.σ_pCO₂,
    yflip = true,
    xlabel = "In-Situ pCO₂ (μatm)",
    ylabel = "Depth (m)",
    label = "pCO₂ ± σ",
    marker = :square,
    linewidth = 2,
    color = :darkred
)

plot(p1, p2, layout = (1, 2), plot_title = "Palaeo CTD Profile Reconstruction")


## =============================================================================
# SECTION 7: Demonstration of stable TA/CO₃⁻² pairs
# ==============================================================================

# Because of the mathematical nature of this pair, when given TA & CO₃⁻² as
# inputs, a single, unique solution is not always given. In many cases, the
# root solvers used for determining pH give two real results. If you are
# concerned that your inputs may produce non-physical or unstable results,
# run this scripts and consult the resultant heat map.

# Testing a range of values including extreme edge cases (μmol/kg)
ta_grid  = range(50.0, 4000.0, length=200)
co3_grid = range(1.0, 800.0, length=200)

# Construct solver for TA and CO3 parameter pair
solve_ta_co3 = CarbonateSystem(:carbon; varying = (:TA, :CO₃), temp_c = 25.0, sal = 35.0)

# Wrapper to catch non-convergence or unphysical pH values
function safe_ph(ta, co3)
    try
        ph = solve_ta_co3(ta, co3).pHtot
        return (isnan(ph) || isinf(ph) || !(0.0 <= ph <= 14.0)) ? NaN : ph
    catch
        return NaN
    end
end

# Broadcast across 2D space (N_ta x N_co3 matrix)
ph_matrix = safe_ph.(ta_grid, co3_grid')

# Render Heatmap
heatmap(
    co3_grid, ta_grid, ph_matrix,
    xlabel = "CO3²⁻ (μmol/kg)",
    ylabel = "Total Alkalinity (μmol/kg)",
    title = "TA vs CO3²⁻ Stability Domain (pH Total)",
    color = :viridis,
    nan_color = :black,
    bg_inside = :black,
    clims = (7.0, 9.5)
)