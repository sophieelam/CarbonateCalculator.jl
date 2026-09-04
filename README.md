# CarbonateCalculator.jl

[![Run Tests and Coverage](https://github.com/sophieelam/CarbonateCalculator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sophieelam/CarbonateCalculator.jl/actions/workflows/ci.yml)

[![codecov](https://codecov.io/github/sophieelam/carbonatecalculator.jl/graph/badge.svg?token=NH92EBMYRY)](https://codecov.io/github/sophieelam/carbonatecalculator.jl)

**CarbonateCalculator.jl** offers a high-performance Julia tool for calculating marine carbonate and boron isotope systems across both modern and paleo ocean settings. Built in reference to [CBSyst](https://github.com/oscarbranson/cbsyst/tree/main), [PyCO2SYS](https://github.com/mvdh7/PyCO2SYS/tree/main), and [CO2SYS-MATLAB](https://github.com/jamesorr/CO2SYS-MATLAB), it optimizes classic oceanographic workflows for Julia:

* **Automatic Differentiation:** Native `ForwardDiff.jl` support enables forward error propagation across all solvers, temperature/pressure corrections, and non-modern seawater compositions (via [`Kgen`](https://github.com/oscarbranson/Kgen)).
* **Palaeo-Compatible:** Powered by pure-Julia equilibrium constants from `Kgen.jl`.
* **Rigorously Validated:** Continuously benchmarked against PyCO2SYS v1.8.3 and CO2SYS-MATLAB v3.2.0.

CarbonateCalculator.jl solves the carbon system from any two parameters, the boron isotope system from any one parameter, or the combined system from flexible parameter sets (see [Equilibrium Constant & Composition Options](#equilibrium-constant--composition-options)). Users can select from 19 carbonate constant formulations, custom boron and calcium composition models, and optional dissociation constants for sulfate, ammonia, and hydrofluoric acid.

---

## Contents

* [Installation](#installation)
* [Supported Input Combinations](#supported-input-combinations)
* [Usage](#usage)
  * [Uncertainty & Error Propagation](#uncertainty--error-propagation)
  * [Recalculating to Collection Conditions](#recalculating-to-collection-conditions-at_collection_conditions)
  * [Vectorization & Batch Processing](#vectorization--batch-processing-broadcasting)
* [Equilibrium Constant & Composition Options](#equilibrium-constant--composition-options)
* [References](#refrences)

---


## Installation

`CarbonateCalculator.jl` relies on [`Kgen`](https://github.com/oscarbranson/Kgen) for palaeo-seawater equilibrium constants. Until registered in the Julia General Registry, install both packages directly via URL:

```julia
using Pkg
Pkg.add(url="[https://github.com/oscarbranson/Kgen](https://github.com/oscarbranson/Kgen)", rev="julia", subdir="julia/Kgen.jl")
Pkg.add(url="[https://github.com/sophieelam/CarbonateCalculator.jl](https://github.com/sophieelam/CarbonateCalculator.jl)")
```

The committed `Manifest.toml` already pins Kgen to that branch, so `Pkg.instantiate()` in a
checkout of this repository is enough on its own. Once Kgen is registered, both steps
collapse into a plain `Pkg.add`.

No Python is required. Constants are calculated by the pure-Julia Kgen, which means
[`ForwardDiff`](https://github.com/JuliaDiff/ForwardDiff.jl) differentiates through them and error propagation works on every method,
including ['MyAMI'](https://github.com/MathisHain/MyAMI). Python is needed only to regenerate the PyCO2SYS comparison baseline in
`test/CO2SYS_tests/`, and is declared in `test/Project.toml` and `test/CondaPkg.toml`
rather than as a dependency of the package.

---

## Supported Input Combinations

| Subsystem Target | Required Inputs | Solved Quantities |
| :--- | :--- | :--- |
| **`:carbon`** | Any 2 of: `TA`, `DIC`, `pH`, `pCO2`, `fCO2`, `CO2`, `HCO3`, `CO3` | Complete carbon speciation, pH, $\Omega_\text{A}$, $\Omega_\text{C}$, $p\text{CO}\_2$ |
| **`:boron`** | Any 1 of: `B_T`, `d11B`, `pH` | Boron speciation and isotopic fractionations ($\delta^{11}\text{B}$) |
| **`:isotopes`** | Any 1 of: `d11B_sw`, `d11B_borate`, `d11B_boric`, `pH` | Boron speciation and isotopic fractionations ($\delta^{11}\text{B}$) |
| **more than one subsystem** | Any valid combination of carbon + boron + boron isotopic parameters | Full carbon and boron system speciation and isotope systematics |

---

## Usage

`CarbonateCalculator` provides two execution workflows tailored to different computational needs:

* **Preset Functions (`carbon_system`, `whole_system`, `boron_system`, `boron_isotopes`):** Best for quick, interactive, single-sample calculations using keyword arguments. Presets handle everyday calculations without requiring solver setup, but do not support Julia vector broadcasting or uncertainty propagation.
* **Compiled Solvers (`CarbonateSystem`):** Built for high-performance batch processing, array broadcasting, and uncertainty propagation. Building a solver object upfront resolves default settings and type parameters once at construction. This enables positional call syntax, seamless vector broadcasting over datasets, and formal uncertainty propagation through `varying_errors`.

#### Preset Functions Example

```julia
using CarbonateCalculator

# Quick single-sample calculation using keyword arguments
res = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 25.0)

println("Calculated pH (total scale): ", res.pHtot)
println("Calculated pCO₂ (μatm):       ", res.pCO₂)
```

#### Compiled Solver Example

```julia
using CarbonateCalculator

# Construct a reusable solver targeting the carbon subsystem
sys = CarbonateSystem(:carbon; varying = (:TA, :DIC))

# Pass parameter values positionally (in the same order as `varying`) to the solver instance
res = sys(2300.0, 2000.0)

println("Calculated pH (total scale): ", res.pHtot)
println("Calculated pCO₂ (μatm):       ", res.pCO₂)
```

### Uncertainty & Error Propagation

`CarbonateCalculator` provides built-in, forward-propagated uncertainty analysis through compiled `CarbonateSystem` solvers. By defining which parameters carry measurement uncertainty via the `varying_errors` argument, the solver automatically calculates and propagates standard deviations to all output quantities.

#### How It Works

1. **Define `varying_errors`:** Specify which inputs have associated uncertainties when constructing your `CarbonateSystem`.
2. **Pass Uncertainties Positionally:** When calling the solver object, provide the standard deviations as additional positional arguments immediately following the parameter values.
3. **Access Results via `.err`:** Output values are accessible as usual, while propagated standard deviations are stored in a dedicated `.err` field on the result structure.

#### Example

```julia
using CarbonateCalculator

# Construct a solver that propagates uncertainties for both TA and DIC
sys = CarbonateSystem(:carbon; 
    varying = (:TA, :DIC), 
    varying_errors = (:TA, :DIC)
)

# Inputs: TA = 2300.0 ± 2.0 μmol/kg, DIC = 2000.0 ± 3.0 μmol/kg
res = sys(2300.0, 2000.0, 2.0, 3.0)

# Access calculated outputs and their propagated uncertainties
println("pH: ", res.pHtot, " ± ", res.err.pHtot)
println("pCO₂: ", res.pCO₂, " ± ", res.err.pCO₂)
```

### Recalculating to Collection Conditions (`at_collection_conditions`)

In marine carbonate chemistry, seawater samples are typically measured in the laboratory at controlled conditions (e.g., 25°C at atmospheric pressure), but researchers need to know the carbonate parameters at the original ocean collection site (e.g., 4°C at 150 bar hydrostatic pressure).

Rather than mixing measurement and collection conditions into a single call with legacy parameters like `T_in` and `T_out`, `CarbonateCalculator` separates this into a clean, two-step process using `at_collection_conditions`.

#### How It Works

1. **Solves the Measurement State:** Calculates the full carbonate system at the laboratory/measurement conditions using supplied inputs (e.g., TA and DIC, or pH and TA).
2. **Preserves Conservative Totals:** Conservative quantities—such as Total Alkalinity (TA), Dissolved Inorganic Carbon (DIC), and Total Boron ($B_T$)—are independent of temperature and pressure. `at_collection_conditions` carries these conserved totals forward.
3. **Re-solves at Collection Conditions:** Re-evaluates equilibrium constants and speciation using the new collection temperature and pressure, returning an updated result set with in-situ pH, saturation states ($\Omega_\text{A}$, $\Omega_\text{C}$) and $p\text{CO}\_2$.

#### Example

```julia
using CarbonateCalculator

# 1. Solve the system at laboratory measurement conditions (25°C, 0 bar)
lab_result = carbon_system(
    TA = 2300.0, 
    DIC = 2000.0, 
    temp_c = 25.0, 
    pres_bar = 0.0
)

# 2. Recalculate parameters at ocean collection conditions (4°C, 150 bar)
field_result = at_collection_conditions(
    lab_result; 
    temp_c = 4.0, 
    pres_bar = 150.0
)

# Compare laboratory vs. in-situ field parameters
println("Lab pH (25°C, 0 bar):       ", lab_result.pHtot)
println("Field pH (4°C, 150 bar):    ", field_result.pHtot)

println("Lab pCO₂ (25°C, 0 bar):     ", lab_result.pCO₂)
println("Field pCO₂ (4°C, 150 bar):  ", field_result.pCO₂)
```

### Vectorization & Batch Processing (Broadcasting)

Julia's native dot broadcasting (`.`) allows compiled `CarbonateSystem` solvers to process large datasets—such as depth profiles, transects, or ocean model outputs—with minimal overhead. Constructing a solver once resolves defaults and type parameters upfront, making array execution fast and memory-efficient.

#### Key Features

* **Zero-Allocation Inner Loops:** Reuses solver metadata across all elements in the input arrays.
* **Flexible Vector Operations:** Combine scalar conditions (e.g., fixed temperature) with vector measurements (e.g., arrays of TA and DIC).
* **Seamless Field Extraction:** Use standard comprehension syntax to extract vectors of calculated output variables.
* **DataFrame Compatibility:** Input and output can be easily integrated with `DataFrames.jl` for structured data analysis.

#### Example

```julia
using CarbonateCalculator

# 1. Construct the solver once
sys = CarbonateSystem(:carbon; varying = (:TA, :DIC))

# 2. Input vectors (e.g., depth profile data)
TA_vec  = [2250.0, 2300.0, 2350.0]
DIC_vec = [1980.0, 2050.0, 2120.0]

# 3. Broadcast the solver across both arrays
results = sys.(TA_vec, DIC_vec)

# 4a. Extract calculated quantities across all samples
pH_profile   = [r.pHtot for r in results]
pCO2_profile = [r.pCO₂ for r in results]

# 4b. Convert output to a DataFrame for structured analysis
using DataFrames
df = DataFrame(results)
profiles = df[:, [:pHtot, :pCO₂]]
```

---

### Equilibrium Constant & Composition Options

Note that by default `K_method="Kgen"`, where *ALL* K values are calculated using [Kgen.jl](https://palaeocarbonatechemistry.github.io/Kgen/), overriding K-specific method flags below.
[**KGen**](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GC011417) provides 'best practice' K values for modern seawater, and provides a polynomial parameterisation for the influence of Mg and Ca variation in palaeo-seawater. 

Equilibrium constant parameterizations and seawater composition models can be specified using `*_method` arguments when constructing a solver or calling preset functions. Full BibTeX entries for all supported methods are available in [`references.bib`](references.bib).

#### Carbonate Dissociation Constants ($K_1, K_2$) — `K_method`

| Key / Method | Reference | Temp (°C) | Salinity | pH Scale & Notes |
| :--- | :--- | :--- | :--- | :--- |
| `"Lueker2000"` | Lueker et al. (2000) | 2 to 35 | 19 to 43 | Total scale; refit of Mehrbach et al. (1973) |
| `"Roy1993"` | Roy et al. (1993) | 0 to 45 | 5 to 45 | Total scale; artificial seawater |
| `"Waters2014"` | Waters et al. (2014) | 0 to 50 | 1 to 50 | Seawater scale; update to Millero et al. (2010) |
| `"SB2020"` | Schockman & Byrne (2020) | 15 to 35 | 19.6 to 40 | Total scale; updated $K_2$ fit |
| `"Millero2010"` | Millero et al. (2010) | 0 to 50 | 1 to 50 | Seawater scale; estuarine & seawater |
| `"Millero2006"` | Millero et al. (2006) | 0 to 50 | 1 to 50 | Seawater scale |
| `"Millero2002"` | Millero et al. (2002) | -1.6 to 35 | 34 to 37 | Seawater scale; field measurements |
| `"MPM2002"` | Mojica Prieto & Millero (2002) | 0 to 45 | 5 to 42 | Seawater scale |
| `"CW2003"` | Cai & Wang (2003) | 2 to 35 | 0 to 49 | Seawater scale (converted from NBS); estuarine |
| `"HM1973"` | Hansson & Mehrbach (1973) | 2 to 35 | 20 to 40 | Seawater scale; refit by Dickson & Millero (1987) |
| `"DM1987"` | Dickson & Millero (1987) | 2 to 35 | 20 to 40 | Seawater scale; refit of Mehrbach et al. (1973) |
| `"Hansson1973"` | Hansson (1973) | 2 to 35 | 20 to 40 | Seawater scale; refit by Dickson & Millero (1987) |
| `"GP1989"` | Goyet & Poisson (1989) | -1 to 40 | 10 to 50 | Seawater scale; artificial seawater |
| `"Mehrbach1973"`| Mehrbach et al. (1973) | 2 to 35 | 19 to 43 | Seawater scale (converted from NBS) |
| `"Papadimitriou2018"` | Papadimitriou et al. (2018)| -6 to 25 | 50 to 100 | Total scale; sea ice & hypersaline waters |
| `"Sulpis2020"` | Sulpis et al. (2020) | 2 to 35 | 19 to 43 | Total scale |
| `"Millero1979"` | Millero (1979) | 0 to 50 | 0 | Pure freshwater formulation ($S = 0$) |

---

#### Total Boron ($B_T$) — `BT_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Uppstrom"` | Uppström (1974) | Standard oceanographic ratio: $B_T = 0.0004157 \cdot S / 35$ |
| `"Lee"` | Lee et al. (2010) | Higher boron formulation: $B_T = 0.0004326 \cdot S / 35$ |
| `"KSK18"` | Kulik et al. (2018) | Non-zero intercept model: $B_T = (10.838 \cdot S + 13.821) \cdot 10^{-6}$ |

---

#### Sulfate & Fluoride Association Constants

**Sulfate Bisulfate ($K_\text{SO4}$) — `KSO4_method`**
| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Dickson"` | Dickson (1990) | Standard oceanographic formulation. |
| `"Kuo"` | Khoo et al. (1977) | Free pH scale formulation. |
| `"WM13"` | Waters & Millero (2013) | Extended temperature and salinity parameterization.

**This setting has no effect if `K_method="Kgen"`.**


**Hydrogen Fluoride ($K_\text{F}$) — `KF_method`**
| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Dickson"` | Dickson & Riley (1979) | Standard oceanographic formulation. |
| `"Perez"` | Pérez & Fraga (1987) | Alternative fit ($S = 10\text{--}40$, $T = 9\text{--}33^\circ\text{C}$).

**This setting has no effect if `K_method="Kgen"`.**


---

#### Minor Species & Calcium Formulations

**Ammonia ($K_\text{NH3}$) — `KNH3_method`**
| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Millero"` | Yao & Millero (1995) | Seawater pH scale formulation. |
| `"Clegg"` | Clegg & Whitfield (1995) | Total pH scale fit ($S = 0\text{--}40$, $T = -2\text{--}40^\circ\text{C}$).

**Calcium Concentration ($[\text{Ca}^{2+}]$) — `Ca_method`**
| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"modern"` | Modern seawater | $[\text{Ca}^{2+}] = 0.01028 \cdot S / 35\ \text{mol/kg-sw}$.
| `"Culkin"` | Culkin (1965) | $[\text{Ca}^{2+}] = 0.01026 \cdot S / 35\ \text{mol/kg-sw}$.
| `"RT67"` | Riley & Tongudai (1967) | $[\text{Ca}^{2+}] = (0.02128 / 40.087) \cdot S / 1.80655\ \text{mol/kg-sw}$.

---

#### Non-Configurable Standard Constants

The following constants use single default parameterizations across the package:

* **$\text{CO}_2$ Solubility ($K_0$):** Weiss (1974)
* **Boric Acid ($K_\text{B}$):** Dickson (1990)
* **Water Dissociation ($K_\text{W}$):** Millero (1995)
* **Phosphoric Acid ($K_\text{P1}, K_\text{P2}, K_\text{P3}$) & Silicic Acid ($K_\text{Si}$):** Yao & Millero (1995)
* **Solubility Products ($K_\text{spA}, K_\text{spC}$):** Mucci (1983)
* **Hydrogen Sulfide ($K_\text{H2S}$):** Millero et al. (1988) / Yao & Millero (1995)
* **Total Fluoride ($F_T$) & Total Sulfate ($S_T$):** Riley (1965) and Morris & Riley (1966)

## References

All equilibrium constant parameterizations, chemical formulations, and isotopic models implemented in `CarbonateCalculator.jl` are compiled in [`references.bib`](references.bib). 

If you use specific constant selections in your work (e.g., Lueker et al., 2000; Lee et al., 2010), you can import `references.bib` directly into Zotero, Mendeley, or BibTeX to cite the original oceanographic literature.