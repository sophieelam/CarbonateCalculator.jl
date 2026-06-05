import PyCO2SYS as pyco
import pandas as pd
import numpy as np
from itertools import product

# --- 1. Aligned Parameters (from your MATLAB script) ---
par_defaults = {1: 2250.0, 2: 2100.0, 3: 8.1, 4: 400.0, 5: 405.0, 6: 200.0, 7: 1800.0, 8: 10.0}

k1k2_opts = list(range(1, 18))
scales_opts = list(range(1, 5))
kso4_opts = [1, 2, 3] # Dickson, Khoo, WM13
kf_opts = [1, 2]
# --- Updated pair_types logic ---
gas_types = {4, 5, 8} # pCO2, fCO2, CO2(aq)
pair_types = [
    (t1, t2) for t1, t2 in product(range(1, 9), range(1, 9)) 
    if t1 != t2 and not (t1 in gas_types and t2 in gas_types)
]
# --- 2. Build Grid ---
grid = list(product(k1k2_opts, scales_opts, kso4_opts, kf_opts, pair_types))
n = len(grid)

k1k2_arr = np.array([g[0] for g in grid])
scale_arr = np.array([g[1] for g in grid])
kso4_arr = np.array([g[2] for g in grid])
kf_arr = np.array([g[3] for g in grid])
p1_type = np.array([g[4][0] for g in grid])
p2_type = np.array([g[4][1] for g in grid])
p1_val = np.array([par_defaults[t] for t in p1_type])
p2_val = np.array([par_defaults[t] for t in p2_type])

# Matched Nutrients
si, phos, nh3, h2s = 10.0, 1.0, 2.0, 3.0

# --- 3. Execute ---
print(f"Calculating {n} rows for PyCO2SYS validation...")
results = pyco.sys(
    par1=p1_val, par2=p2_val, par1_type=p1_type, par2_type=p2_type,
    salinity=33.1, temperature=24.0, pressure=1.0,
    total_silicate=si, total_phosphate=phos, total_ammonia=nh3, total_sulfide=h2s,
    opt_k_carbonic=k1k2_arr, opt_pH_scale=scale_arr,
    opt_k_bisulfate=kso4_arr, opt_k_fluoride=kf_arr, opt_total_borate=1
)

# --- 4. Export ---
df = pd.DataFrame({
    'PAR1': p1_val, 'PAR2': p2_val, 'PAR1TYPE': p1_type, 'PAR2TYPE': p2_type,
    'K1K2CONSTANTS': k1k2_arr, 'pHSCALEIN': scale_arr, 'KSO4CONSTANT': kso4_arr,
    'KFCONSTANT': kf_arr, 'BORON': 1, 'TEMPIN': 24.0, 'SAL': 33.1, 'PRESIN': 10.0,
    'PO4': phos, 'SI': si, 'H2S': h2s, 'NH4': nh3,
    'DIC_out': results['dic'], 'pH_out': results['pH'], 
    'fCO2_out': results['fCO2'], 'CO3_out': results['CO3']
})
df.to_csv("results/compare_PyCO2SYS_v1_8_3.csv", index=False)
print("Digital twin baseline complete!")



import PyCO2SYS as pyco2

# Row 18 inputs
results = pyco2.sys(
    par1=8.1, par1_type=3,    # pH
    par2=405.0, par2_type=5,  # fCO2
    temperature=24.0,
    salinity=33.1,
    pressure=10.0,            # 10 dbar (equivalent to your 1.0 bar gauge)
    opt_k_carbonic=1,         # Roy 1993
    opt_pH_scale=1,           # Total Scale
    opt_k_fluoride=2,         # Perez (from your PyCO2SYS mapping)
    opt_k_sulfate=1,          # Dickson
    opt_b_value=2             # Lee (from your PyCO2SYS mapping)
)

print(f"PyCO2SYS Target K0: {results['k_CO2']}")
print(f"PyCO2SYS Target K1: {results['k_carbonic_1']}")
print(f"PyCO2SYS Target K2: {results['k_carbonic_2']}")