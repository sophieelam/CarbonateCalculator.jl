"""Generate the PyCO2SYS v1.8.3 baselines for test/external/compare_*.jl.

    pip install "PyCO2SYS==1.8.3.4"
    python test/external/generate_pyco2sys_baseline.py

Writes two files to test/external/results/:

* compare_PyCO2SYS_v1_8_3.csv: every K1/K2, pH-scale, KSO4 and KF option and
  every valid input pair, at one condition (S 33.1, T 24 °C, 10 dbar), with output
  conditions 12 °C and 1647 dbar (the same conditions as the MATLAB baseline).
* compare_PyCO2SYS_v1_8_3_conditions.csv: a sweep over temperature, salinity and pressure
  (0 to 11 000 dbar) for three K1/K2 options of each native pH scale, on the total and
  free scales. This is the file that tests the pressure corrections.

Columns use CO2SYS-MATLAB names (TAlk, TCO2, pHin, pCO2in, ..., K1input, ...), so one
reader in agreement_common.jl serves both references. PRESIN and PRESOUT are in dbar, as
in CO2SYS. Constants are on the input pH scale, except KS and KF (free scale), as in
CO2SYS; KSPCinput and KSPAinput (solubility products) are extra columns CO2SYS lacks.

Input-type codes in PAR1TYPE/PAR2TYPE follow PyCO2SYS: 6 = CO3, 7 = HCO3 (CO2SYS v3
swaps these two).
"""

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import PyCO2SYS as pyco

assert pyco.__version__.startswith("1.8.3"), f"PyCO2SYS {pyco.__version__}: need 1.8.3"

RESULTS = Path(__file__).parent / "results"
PAR_DEFAULTS = {1: 2250.0, 2: 2100.0, 3: 8.1, 4: 400.0, 5: 405.0, 6: 200.0, 7: 1800.0, 8: 10.0}
GAS = {4, 5, 8}  # pCO2, fCO2, CO2(aq): no two of these make a valid pair
PAIRS = [(a, b) for a, b in product(range(1, 9), range(1, 9))
         if a != b and not (a in GAS and b in GAS)]
NUTRIENTS = dict(SI=10.0, PO4=1.0, NH4=2.0, H2S=3.0)

# PyCO2SYS result key -> CO2SYS-MATLAB column (input conditions; "_out" keys -> "...out")
OUTPUTS = {
    "alkalinity": "TAlk", "dic": "TCO2", "pH": "pHin", "pCO2": "pCO2in", "fCO2": "fCO2in",
    "HCO3": "HCO3in", "CO3": "CO3in", "aqueous_CO2": "CO2in",
    "saturation_calcite": "OmegaCAin", "saturation_aragonite": "OmegaARin",
    "revelle_factor": "RFin",
    "pH_total": "pHinTOTAL", "pH_sws": "pHinSWS", "pH_free": "pHinFREE", "pH_nbs": "pHinNBS",
}
OUTPUTS_OUT = {
    "pH_out": "pHout", "pCO2_out": "pCO2out", "fCO2_out": "fCO2out", "HCO3_out": "HCO3out",
    "CO3_out": "CO3out", "aqueous_CO2_out": "CO2out",
    "saturation_calcite_out": "OmegaCAout", "saturation_aragonite_out": "OmegaARout",
    "revelle_factor_out": "RFout",
    "pH_total_out": "pHoutTOTAL", "pH_sws_out": "pHoutSWS", "pH_free_out": "pHoutFREE",
    "pH_nbs_out": "pHoutNBS",
}
CONSTANTS = {
    "k_CO2": "K0input", "k_carbonic_1": "K1input", "k_carbonic_2": "K2input",
    "k_water": "KWinput", "k_borate": "KBinput", "k_fluoride": "KFinput",
    "k_bisulfate": "KSinput", "k_phosphoric_1": "KP1input", "k_phosphoric_2": "KP2input",
    "k_phosphoric_3": "KP3input", "k_silicate": "KSiinput", "k_sulfide": "KH2Sinput",
    "k_ammonia": "KNH3input", "k_calcite": "KSPCinput", "k_aragonite": "KSPAinput",
}


def run(grid: pd.DataFrame) -> pd.DataFrame:
    p1 = grid.PAR1TYPE.map(PAR_DEFAULTS).to_numpy()
    p2 = grid.PAR2TYPE.map(PAR_DEFAULTS).to_numpy()
    # No nutrients with the freshwater constants (option 8): CO2SYS zeroes their equilibrium
    # constants there, which CarbonateCalculator.jl does not (a documented difference), so
    # nutrients would hide the freshwater carbonate chemistry this row is meant to test.
    fresh = grid.K1K2CONSTANTS.to_numpy() == 8
    nut = {k: np.where(fresh, 0.0, v) for k, v in NUTRIENTS.items()}
    r = pyco.sys(
        par1=p1, par2=p2, par1_type=grid.PAR1TYPE.to_numpy(), par2_type=grid.PAR2TYPE.to_numpy(),
        salinity=grid.SAL.to_numpy(), temperature=grid.TEMPIN.to_numpy(),
        pressure=grid.PRESIN.to_numpy(), temperature_out=grid.TEMPOUT.to_numpy(),
        pressure_out=grid.PRESOUT.to_numpy(),
        total_silicate=nut["SI"], total_phosphate=nut["PO4"],
        total_ammonia=nut["NH4"], total_sulfide=nut["H2S"],
        opt_k_carbonic=grid.K1K2CONSTANTS.to_numpy(), opt_pH_scale=grid.pHSCALEIN.to_numpy(),
        opt_k_bisulfate=grid.KSO4CONSTANT.to_numpy(), opt_k_fluoride=grid.KFCONSTANT.to_numpy(),
        opt_total_borate=grid.BORON.to_numpy(),
    )
    out = grid.copy()
    out["PAR1"], out["PAR2"] = p1, p2
    for k, v in nut.items():
        out[k] = v
    for table in (OUTPUTS, OUTPUTS_OUT, CONSTANTS):
        for key, col in table.items():
            out[col] = np.broadcast_to(np.asarray(r[key], dtype=float), len(out))
    return out


def options_sweep() -> pd.DataFrame:
    # Boron option 1 only, as before: the MATLAB baseline covers both boron options.
    rows = product(range(1, 18), range(1, 5), (1, 2, 3), (1, 2), (1,), PAIRS)
    g = pd.DataFrame(rows, columns=["K1K2CONSTANTS", "pHSCALEIN", "KSO4CONSTANT",
                                    "KFCONSTANT", "BORON", "PAIR"])
    # Freshwater K1/K2 (option 8) only at S = 0; everything else at S = 33.1
    g["SAL"] = np.where(g.K1K2CONSTANTS == 8, 0.0, 33.1)
    g["TEMPIN"], g["PRESIN"], g["TEMPOUT"], g["PRESOUT"] = 24.0, 10.0, 12.0, 1647.0
    return g


def conditions_sweep() -> pd.DataFrame:
    # One total-scale, one seawater-scale and one more total-scale K1/K2 fit, so that every
    # route through the pressure correction is exercised: Roy 1993 (total), Lueker 2000
    # (total), Millero 2010 (seawater).
    rows = product((0.0, 10.0, 25.0, 35.0), (20.0, 33.1, 38.0),
                   (0.0, 1000.0, 5000.0, 11000.0), (1, 10, 14), (1, 3), PAIRS)
    g = pd.DataFrame(rows, columns=["TEMPIN", "SAL", "PRESIN", "K1K2CONSTANTS",
                                    "pHSCALEIN", "PAIR"])
    g["KSO4CONSTANT"], g["KFCONSTANT"], g["BORON"] = 1, 1, 2
    g["TEMPOUT"], g["PRESOUT"] = 2.0, 4000.0
    return g


def finish(g: pd.DataFrame) -> pd.DataFrame:
    g["PAR1TYPE"] = [p[0] for p in g.PAIR]
    g["PAR2TYPE"] = [p[1] for p in g.PAIR]
    return run(g.drop(columns="PAIR"))


if __name__ == "__main__":
    RESULTS.mkdir(exist_ok=True)
    for name, grid in [("compare_PyCO2SYS_v1_8_3.csv", options_sweep()),
                       ("compare_PyCO2SYS_v1_8_3_conditions.csv", conditions_sweep())]:
        df = finish(grid)
        df.to_csv(RESULTS / name, index=False, float_format="%.17g")
        print(f"{name}: {len(df):,} rows")
