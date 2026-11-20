# CarbonateCalculator.jl

[![Run Tests and Coverage](https://github.com/sophieelam/CarbonateCalculator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sophieelam/CarbonateCalculator.jl/actions/workflows/ci.yml)

[![codecov](https://codecov.io/github/sophieelam/carbonatecalculator.jl/graph/badge.svg?token=NH92EBMYRY)](https://codecov.io/github/sophieelam/carbonatecalculator.jl)

**CarbonateCalculator.jl** offers a high-performance Julia tool for calculating marine carbonate and boron isotope systems across both modern and paleo ocean settings. Built in reference to [CBSyst](https://github.com/oscarbranson/cbsyst/tree/main), [PyCO2SYS](https://github.com/mvdh7/PyCO2SYS/tree/main), and [CO2SYS-MATLAB](https://github.com/jamesorr/CO2SYS-MATLAB), it optimizes classic oceanographic workflows for Julia:

* **Automatic Differentiation:** Native `ForwardDiff.jl` support enables forward error propagation across all solvers, temperature/pressure corrections, and non-modern seawater compositions (via [`Kgen`](https://github.com/oscarbranson/Kgen)).
* **Palaeo-Compatible:** Powered by pure-Julia equilibrium constants from `Kgen.jl`.
* **Rigorously Validated:** Continuously benchmarked against PyCO2SYS v1.8.3 and CO2SYS-MATLAB v3.2.0.

CarbonateCalculator.jl solves the carbon system from any two parameters, the boron isotope system from any one parameter, or the combined system from flexible parameter sets (see [Chemical Constants](#chemical-constants)). Users can select from 19 carbonate constant formulations, custom boron and calcium composition models, and optional dissociation constants for sulfate, ammonia, and hydrofluoric acid.

---

## Contents

* [Installation](#installation)
* [Supported Input Combinations](#supported-input-combinations)
* [Usage](#usage)
  * [Uncertainty & Error Propagation](#uncertainty--error-propagation)
  * [Changing Temperature & Pressure](#changing-temperature--pressure-adjust)
  * [Input Pairs With Two Solutions](#input-pairs-with-two-solutions)
  * [Vectorization & Batch Processing](#vectorization--batch-processing-broadcasting)
* [Reference: Arguments & Results](#reference-arguments--results)
* [Chemical Constants](#chemical-constants)
  * [Choosing the constants (`K_settings`)](#choosing-the-constants-k_settings)
  * [Warnings](#warnings)
  * [Moving from `K_method`](#moving-from-k_method)
* [References](#references)

---


## Installation

`CarbonateCalculator.jl` relies on [`Kgen`](https://github.com/oscarbranson/Kgen) for palaeo-seawater equilibrium constants. Until registered in the Julia General Registry, install both packages directly via URL:

```julia
using Pkg
Pkg.add(url="https://github.com/oscarbranson/Kgen", rev="julia", subdir="julia/Kgen.jl")
Pkg.add(url="https://github.com/sophieelam/CarbonateCalculator.jl")
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
| **`:carbon`** | Any 2 of: `TA`, `DIC`, pH (`pHtot`/`pHsws`/`pHfree`/`pHNBS`), `CO₂` (or `pCO₂`/`fCO₂`), `HCO₃`, `CO₃` (or `ΩA`/`ΩC`) | Complete carbon speciation, pH, $\Omega_\text{A}$, $\Omega_\text{C}$, $p\text{CO}_\text{2}$ |
| **`:boron`** | Any 1 of: pH, `BOH₃`, `BOH₄` (total boron `BT` comes from salinity unless given) | Boron speciation |
| **`:isotopes`** | Any 1 of: pH, `δBOH₄` (or `ABOH₄`), `δBOH₃` (or `ABOH₃`) (`δBT` is modern seawater unless given) | Boron isotope speciation ($\delta^{11}\text{B}$) |
| **more than one subsystem** | Any valid combination of carbon + boron + boron isotopic parameters | Full carbon and boron system speciation and isotope systematics |

Two of the carbon pairs, `TA` with `CO₃` and `DIC` with `HCO₃`, usually have **two valid solutions**; see [Input Pairs With Two Solutions](#input-pairs-with-two-solutions). Every argument and result is listed, with units, in [Reference: Arguments & Results](#reference-arguments--results).

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

### Changing Temperature & Pressure (`adjust`)

`adjust` re-solves a result at a different temperature and/or pressure: the same water, with the same TA, DIC, totals and composition, at new conditions. A common use is reporting samples measured in the laboratory (e.g. 25 °C at atmospheric pressure) at the conditions they were collected at (e.g. 4 °C at 150 bar), but it applies equally to normalising a dataset to one temperature, or following a water mass as it warms or sinks. The name follows PyCO2SYS v2, whose equivalent is also `adjust`.

Rather than mixing two sets of conditions into one call with legacy parameters like `T_in` and `T_out`, a calculation describes one set of conditions, and `adjust` moves its result to another.

#### How It Works

1. **Start from a solved result.** `adjust` takes the result of a calculation at the original conditions (from `carbon_system`, `whole_system`, or a `CarbonateSystem`), which already holds TA and DIC.
2. **Carry the conservative quantities.** TA, DIC and total boron ($B_T$) do not change with temperature or pressure, so they are carried across unchanged, along with the seawater composition and every method choice.
3. **Solve once at the new conditions.** The equilibrium constants are recalculated at the new temperature and pressure, and the system is solved from TA and DIC, giving the new pH, saturation states ($\Omega_\text{A}$, $\Omega_\text{C}$) and $p\text{CO}_\text{2}$.

The original conditions are not solved again: TA and DIC come straight from the result. The one exception is uncertainty propagation when the original pair was *not* TA and DIC (pH and TA, say). There, TA or DIC was derived using constants at the original conditions, so that solve is repeated inside the differentiation to carry its uncertainty through. When the original pair is TA and DIC, it isn't needed and is skipped.

> **Tip:** if you measured TA and DIC and only need the system at one set of conditions, you don't need `adjust` at all; solve directly there: `carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 4.0, pres_bar = 150.0)`.

#### Example

```julia
using CarbonateCalculator

# 1. Solve the system at laboratory measurement conditions (25°C, 0 bar)
lab_result = carbon_system(
    TA = 2300.0, 
    pHtot = 7.9,
    temp_c = 25.0, 
    pres_bar = 0.0
)

# 2. Re-solve the same water at ocean collection conditions (4°C, 150 bar)
field_result = adjust(
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

`adjust` also broadcasts over many results (`adjust.(results, temps, pressures)`), and takes `σ_temp_c` and `σ_pres_bar` for uncertainty in the new conditions. `at_collection_conditions`, its former name, still works but is deprecated.

### Input Pairs With Two Solutions

Two input pairs do not always determine a unique system, and it is possible to have **two pH solutions that are both physically valid** ([Humphreys et al., 2022](https://gmd.copernicus.org/articles/15/15/2022/), Sect. 5.2 and Fig. 7).

| Pair | Setting | Default | The other root |
| :--- | :--- | :--- | :--- |
| `TA` with `CO₃` (or `ΩA`/`ΩC`) | `CO3_root` | `"low"`: the lower-pH, seawater-like solution | `"high"`: a hydroxide-dominated solution with little DIC |
| `DIC` with `HCO₃` | `HCO3_root` | `"high"`: the higher-pH solution | `"low"`: the CO₂-rich solution on the other side of the bicarbonate peak |

The defaults are the conventional choices and match PyCO2SYS. For example, the case in Humphreys et al. (2022, Fig. 7):

```julia
carbon_system(TA = 2300.0, CO₃ = 120.0, temp_c = 15.0).pHtot                     # ≈ 7.9  (DIC ≈ 2144)
carbon_system(TA = 2300.0, CO₃ = 120.0, temp_c = 15.0, CO3_root = "high").pHtot  # ≈ 10.8 (DIC ≈ 122)
```

A TA that is too low for its `CO₃` has no solution at all and raises an error; a DIC–`HCO₃` pair with no solution returns `NaN`. Every other pair has exactly one solution and ignores both settings.

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

## Reference: Arguments & Results

### Arguments

Every function that solves a system (`carbon_system`, `whole_system`, `boron_system`, `boron_isotopes`, and `CarbonateSystem`) takes these names. Each scope accepts the ones relevant to it; an argument outside the scope is an error naming the scope that would accept it. Concentrations are in the unit set by `unit` (µmol/kg by default) unless marked otherwise.

| Argument | Meaning | Unit | Default | Scope |
| :--- | :--- | :--- | :--- | :--- |
| **Carbonate system** (give two) | | | | |
| `TA` | Total alkalinity | `unit` | — | carbon |
| `DIC` | Dissolved inorganic carbon | `unit` | — | carbon |
| `pHtot`, `pHsws`, `pHfree`, `pHNBS` | pH on the total, seawater, free or NBS scale (one measurement; give one) | — | — | all |
| `CO₂` | Aqueous CO₂ | `unit` | — | carbon |
| `pCO₂`, `fCO₂` | CO₂ partial pressure, fugacity (stand in for `CO₂`) | µatm | — | carbon |
| `HCO₃` | Bicarbonate | `unit` | — | carbon |
| `CO₃` | Carbonate ion | `unit` | — | carbon |
| `ΩA`, `ΩC` | Saturation state, aragonite or calcite (stand in for `CO₃`) | — | — | carbon |
| **Boron** | | | | |
| `BOH₃`, `BOH₄` | Boric acid, borate | `unit` | — | boron |
| **Boron isotopes** | | | | |
| `δBOH₃`, `δBOH₄` | δ¹¹B of boric acid, borate | ‰ | — | isotopes |
| `ABOH₃`, `ABOH₄` | ¹¹B fractional abundance of boric acid, borate (alternatives to δ) | — | — | isotopes |
| `δBT` | δ¹¹B of total boron (seawater) | ‰ | modern seawater, 39.61 | isotopes |
| `ABT` | ¹¹B fractional abundance of total boron | — | from `δBT` | isotopes |
| `alphaB` | Boron isotope fractionation factor, α(B(OH)₃–B(OH)₄⁻) | — | Klochko et al. (2006) | isotopes |
| **Conditions** | | | | |
| `temp_c` | Temperature | °C | 25.0 | all |
| `sal` | Practical salinity | — | 35.0 | all |
| `pres_bar` | Hydrostatic pressure (0 at the surface) | bar | 0.0 | all |
| **Totals & composition** | | | | |
| `BT` | Total boron | `unit` | from salinity (`K_settings.BT_method`) | all |
| `ST` | Total sulphate | mol/kg | from salinity | all |
| `FT` | Total fluoride | mol/kg | from salinity | all |
| `PT`, `SiT` | Total phosphate, silicate | `unit` | 0.0 | carbon |
| `H2ST`, `NH4T` | Total hydrogen sulphide, ammonia | `unit` | 0.0 | carbon |
| `Ca`, `Mg` | Calcium, magnesium, given at S = 35 | mol/kg | modern seawater (`K_settings.Ca_method`) | all |
| **Settings** | | | | |
| `unit` | Concentration unit: `"mol"`, `"mmol"`, `"umol"`, `"nmol"`, `"pmol"`, `"fmol"` (per kg) | — | `"umol"` | all |
| `K_settings` | Which parameterisation to use for each constant, as a NamedTuple naming only what differs from the defaults, e.g. `(K1_K2_method = "Roy 1993", KF_method = "Perez")`, or a `KSettings` | — | `KSettings()` — see [Choosing the constants](#choosing-the-constants-k_settings) | all |
| `CO3_root` | Which solution for TA with `CO₃`/`ΩA`/`ΩC`: `"low"` or `"high"` pH | — | `"low"` | all |
| `HCO3_root` | Which solution for DIC with `HCO₃`: `"low"` or `"high"` pH | — | `"high"` | all |
| `Ks` | A pre-computed constant bundle from `calculate_constants`, used instead of computing one | — | — | all |

`CarbonateSystem` additionally takes `varying` (the parameters given positionally) and `varying_errors` (those with uncertainties, given as σ in the same unit as the parameter). `adjust` takes `temp_c`, `pres_bar`, `σ_temp_c` and `σ_pres_bar`.

### Results

A result's values are read as fields, e.g. `res.pHtot`, and uncertainties as `res.err.pHtot` (or the `σ_pHtot` column of a `DataFrame`). Which keys are present depends on the scope.

| Key | Meaning | Unit | Scope |
| :--- | :--- | :--- | :--- |
| `pHtot`, `pHsws`, `pHfree`, `pHNBS` | pH on each scale | — | all |
| `temp_c`, `sal`, `pres_bar` | Conditions | °C, —, bar | all |
| `BT` | Total boron | `unit` | all |
| `ST`, `FT` | Total sulphate, fluoride | mol/kg | all |
| `Ca`, `Mg` | Calcium, magnesium, at S = 35 | mol/kg | all |
| `fH` | Hydrogen ion activity coefficient, for the NBS scale | — | all |
| `TA`, `DIC` | Total alkalinity, dissolved inorganic carbon | `unit` | carbon |
| `CO₂`, `HCO₃`, `CO₃` | Carbonate species | `unit` | carbon |
| `pCO₂`, `fCO₂` | CO₂ partial pressure, fugacity | µatm | carbon |
| `ΩA`, `ΩC` | Saturation state, aragonite and calcite | — | carbon |
| `H` | [H⁺] on the total scale | mol/kg | carbon |
| `Hfree` | Free [H⁺] | `unit` | carbon |
| `CAlk`, `BAlk`, `OH`, `PAlk`, `SiAlk`, `Alk_NH3`, `Alk_H2S` | Alkalinity contributions: carbonate, borate, hydroxide, phosphate, silicate, ammonia, sulphide | `unit` | carbon |
| `HSO₄`, `HF` | Bisulphate, hydrogen fluoride (negative alkalinity terms) | `unit` | carbon |
| `PT`, `SiT`, `H2ST`, `NH4T` | Nutrient totals, as used | `unit` | carbon |
| `BOH₃`, `BOH₄` | Boric acid, borate | `unit` | boron |
| `δBT`, `δBOH₃`, `δBOH₄` | δ¹¹B of total boron, boric acid, borate | ‰ | isotopes |
| `ABT`, `ABOH₃`, `ABOH₄` | ¹¹B fractional abundances | — | isotopes |
| `alphaB` | Fractionation factor used | — | isotopes |

Besides these values, a result also carries `res.Ks` (the equilibrium constants used, in mol/kg on the total pH scale, with `KS` and `KF` on the free scale), `res.unit`, and `res.settings` (`K_settings`, `unit`, `CO3_root` and `HCO3_root`; `res.settings.K_settings.K1_K2_method` is the K₁/K₂ parameterisation used). `revelle_factor(res)` and `calc_gradient` give buffer factors.

## Chemical Constants

All constants and composition functions are specified in [`constants.jl`](src/constants.jl).

### Seawater Composition

For conservative elements, composition is estimated from salinity using the following formulations:

#### Total Boron ($B_T$) — `BT_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Uppstrom"`* | Uppström (1974) | Standard oceanographic ratio: $B_T = 0.0004157 \cdot S / 35$ |
| `"Lee"` | Lee et al. (2010) | Higher boron formulation: $B_T = 0.0004326 \cdot S / 35$ |
| `"KSK18"` | Kulik et al. (2018) | Non-zero intercept model: $B_T = (10.838 \cdot S + 13.821) \cdot 10^{-6}$ |
| `"Culkin"` | Culkin (1965), as used by GEOSECS and Peng et al. | $B_T = 0.0004106 \cdot S / 35$, about 1 % below Uppström |
| `"freshwater"` | — | $B_T = 0$. Intended for `K1_K2_method = "Millero 1979"`; warns with anything else |

#### Calcium Concentration ($[\text{Ca}^{2+}]$) — `Ca`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"modern"`* | Modern seawater from `KGen` | $[\text{Ca}^{2+}] = 0.01028 \cdot S / 35\ \text{mol/kg-sw}$.
| `"Culkin"` | Culkin (1965) | $[\text{Ca}^{2+}] = 0.01026 \cdot S / 35\ \text{mol/kg-sw}$.
| `"RT67"` | Riley & Tongudai (1967) | $[\text{Ca}^{2+}] = (0.02128 / 40.087) \cdot S / 1.80655\ \text{mol/kg-sw}$.
| `"[a number]"` | User-settable | The concentration of calcium in the seawater in $\text{mol/kg-sw}$ at 35.0 salinity. |

#### Magnesium Concentration ($[\text{Mg}^{2+}]$) — `Mg`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"modern"`* | Modern seawater from `KGen` | $[\text{Mg}^{2+}] = 0.0128 \cdot S / 35\ \text{mol/kg-sw}$.
| `"[a number]"` | User-settable | The concentration of magnesium in the seawater in $\text{mol/kg-sw}$ at 35.0 salinity. |

#### Total Sulphate ($S_T$)

Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966

#### Total Fluoride ($F_T$)

Riley, J. P., Deep-Sea Research 12:219-220, 1965

### Equilibrium Constants

Every parameterisation choice is made through `K_settings` (see [Choosing the constants](#choosing-the-constants-k_settings) below); the field each table belongs to is in its heading, and `*` marks the default. Full BibTeX entries for all supported methods are available in [`references.bib`](references.bib).

[**KGen**](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GC011417), the default `K1_K2_method`, provides 'best practice' K values for modern seawater, and a polynomial parameterisation for the influence of Mg and Ca variation in palaeo-seawater. Under KGen, Kgen.jl supplies K₀, K₁, K₂, K_B, K_W, K_SO4, K_F, the phosphate and silicate constants and both solubility products, already pressure-corrected, so `KSO4_method` and `KF_method` have no effect and `pressure_method` reaches only K_H2S and K_NH3. Choosing them alongside KGen warns.

#### Carbonate ($K_1, K_2$) — `K1_K2_method`

Each parameterisation was fitted over the temperature and salinity ranges below. A calculation outside them still runs, but warns (see [Warnings](#warnings)). These ranges are the registry in [`constants.jl`](src/constants.jl) (`SPECIATION.K1_K2_method`), which is what the warning reads; ranges are from the CO2SYS v3.2 header unless noted.

| Key | Reference | Temp (°C) | Salinity | Native pH scale & notes | Range source |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `"KGen"`* | Whiteford et al. (2025) | −2 to 50 | 0 to 50 | Total; modern & palaeo seawater | this README (to be confirmed — Kgen's own documentation gives none) |
| `"Lueker 2000"` | Lueker et al. (2000) | 2 to 35 | 19 to 43 | Total; refit of Mehrbach et al. (1973) | CO2SYS option 10 |
| `"Roy 1993"` | Roy et al. (1993) | 0 to 45 | 5 to 45 | Total; artificial seawater | CO2SYS option 1 |
| `"GP 1989"` | Goyet & Poisson (1989) | −1 to 40 | 10 to 50 | Seawater; artificial seawater | CO2SYS option 2 |
| `"Hansson 1973"` | Hansson (1973), refit by Dickson & Millero (1987) | 2 to 35 | 20 to 40 | Seawater | CO2SYS option 3 |
| `"DM 1987"` | Mehrbach et al. (1973), refit by Dickson & Millero (1987) | 2 to 35 | 20 to 40 | Seawater | CO2SYS option 4 |
| `"HM 1973"` | Hansson and Mehrbach combined, Dickson & Millero (1987) | 2 to 35 | 20 to 40 | Seawater | CO2SYS option 5 |
| `"Mehrbach 1973"` | Mehrbach et al. (1973) | 2 to 35 | 19 to 43 | NBS, converted with `fH`; has companion sets — see below | CO2SYS options 6, 7 |
| `"Millero 1979"` | Millero (1979) | 0 to 50 | 0 only | NBS; pure water; has a companion set — see below | CO2SYS option 8 |
| `"CW 2003"` | Cai & Wang (1998) | 2 to 35 | 0 to 49 | NBS, converted to seawater scale; estuarine | CO2SYS option 9 |
| `"MPM 2002"` | Mojica Prieto & Millero (2002) | 0 to 45 | 5 to 42 | Seawater | CO2SYS option 11 |
| `"Millero 2002"` | Millero et al. (2002) | −1.6 to 35 | 34 to 37 | Seawater; field measurements | CO2SYS option 12 |
| `"Millero 2006"` | Millero et al. (2006) | 0 to 50 | 1 to 50 | Seawater | CO2SYS option 13 |
| `"Millero 2010"` | Millero (2010) | 0 to 50 | 1 to 50 | Seawater; estuarine & seawater | CO2SYS option 14 |
| `"Waters 2014"` | Waters et al. (2014) | 0 to 50 | 1 to 50 | Seawater; update to Millero (2010) | CO2SYS option 15 |
| `"Sulpis 2020"` | Sulpis et al. (2020) | −1.7 to 32 | 31 to 38 | Total; field data | CO2SYS option 16 |
| `"SB 2020"` | Schockman & Byrne (2021) | 15 to 35 | 19 to 41 | Total; updated $K_2$ fit | CO2SYS option 17 |
| `"Papadimitriou 2018"` | Papadimitriou et al. (2018) | −6 to 25 | 33 to 100 | Total; sea ice & hypersaline brines | the paper |

Outside the union of all these ranges (−6 to 50 °C, salinity 0 to 100) no parameterisation applies, and every calculation warns regardless of the one selected.

#### Companion sets

Two K₁/K₂ parameterisations were published, or are conventionally used, with particular choices for other constants. Selecting one of them does **not** change anything else; instead, if your `K_settings` does not match one of its sets, you get a warning naming the sets and what you currently have selected.

| `K1_K2_method` | Set | `K_settings` |
| :--- | :--- | :--- |
| `"Mehrbach 1973"` | GEOSECS (Takahashi et al., 1982; CO2SYS option 6) | `(K1_K2_method = "Mehrbach 1973", fH_method = "Takahashi 1982", BT_method = "Culkin", Ca_method = "Culkin")` |
| `"Mehrbach 1973"` | Peng et al. (1987; CO2SYS option 7) | `(K1_K2_method = "Mehrbach 1973", fH_method = "Peng 1987", BT_method = "Culkin", Ca_method = "Culkin")` |
| `"Millero 1979"` | Freshwater (CO2SYS option 8) | `(K1_K2_method = "Millero 1979", fH_method = "freshwater", BT_method = "freshwater", pressure_method = "Millero 1983")` |

(`fH_method = "Takahashi 1982"` is the default, so the GEOSECS set can leave it out.) The companion-only options — `fH_method = "Peng 1987"`, `fH_method = "freshwater"`, `BT_method = "freshwater"` and `pressure_method = "Millero 1983"` — warn in the other direction when used with any other K₁/K₂ parameterisation.

CO2SYS changes more than this under options 6–8, and those parts are **not implemented** here, so results under these sets agree with CO2SYS for the carbonate system but not in every detail:

* GEOSECS and Peng et al. (options 6, 7): K_B from Li et al. (1969) on the NBS scale; K_W of zero (GEOSECS) or Millero (1979)'s seawater K_W (Peng); zero phosphate and silicate constants (GEOSECS) or Peng's fixed values; the GEOSECS pressure corrections and solubility products; fCO₂ equal to pCO₂; and Peng's convention of including phosphate in alkalinity.
* Freshwater (option 8): CO2SYS also forces salinity to 0. Here salinity is whatever you give, and `"Millero 1979"` warns if it is not 0. CO2SYS also sets K_B, the phosphate and silicate constants, K_H2S and K_NH4 to zero, so phosphate, silicate, sulfide and ammonia add no alkalinity; here they are computed as usual. With no boron or nutrients the two agree; with `PT`, `SiT`, `H2ST` or `NH4T` given in fresh water, this package includes their alkalinity and CO2SYS does not.

#### Sulfate ($K_\text{SO4}$, $K_\text{HSO4}$) — `KSO4_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Dickson"`* | Dickson (1990) | Standard oceanographic formulation. |
| `"Khoo"` | Khoo et al. (1977) | Free pH scale formulation. |
| `"WM13"` | Waters & Millero (2013) | Extended temperature and salinity parameterization.

Has no effect under `K1_K2_method = "KGen"`, and warns if chosen with it.

#### Hydrogen Fluoride ($K_\text{F}$) — `KF_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Dickson"`* | Dickson & Riley (1979) | Standard oceanographic formulation. |
| `"Perez"` | Pérez & Fraga (1987) | Alternative fit ($S = 10\text{--}40$, $T = 9\text{--}33^\circ\text{C}$).

Has no effect under `K1_K2_method = "KGen"`, and warns if chosen with it.

#### Ammonia ($K_\text{NH3}$) — `KNH3_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Millero"`* | Yao & Millero (1995) | Seawater pH scale formulation. |
| `"Clegg"` | Clegg & Whitfield (1995) | Total pH scale fit ($S = 0\text{--}40$, $T = -2\text{--}40^\circ\text{C}$).

#### NBS activity coefficient ($f_\text{H}$) — `fH_method`

$f_\text{H}$ relates the NBS pH scale to the seawater scale. It converts `pHNBS` and the NBS-fitted K₁/K₂ (Mehrbach 1973, Millero 1979, Cai & Wang).

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Takahashi 1982"`* | Takahashi et al. (1982), GEOSECS Pacific Expedition v. 3 | The CO2SYS default for every option but 7 and 8. |
| `"Peng 1987"` | Peng et al. (1987) | Intended for `"Mehrbach 1973"` (the Peng et al. set). |
| `"freshwater"` | — | $f_\text{H} = 1$. Intended for `"Millero 1979"`. |

#### Pressure corrections — `pressure_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Millero 1995"`* | Millero (1995) | Partial molal volumes and compressibilities for every constant. |
| `"Millero 1983"` | Millero (1983) | Freshwater coefficients for K₁, K₂ and K_W (K_B uncorrected), the rest as Millero (1995). Intended for `"Millero 1979"`. |

Under KGen this reaches only K_H2S and K_NH3; Kgen corrects its own constants.

Both coefficient sets are for constants on the seawater pH scale, so every correction is applied there, as in CO2SYS and PyCO2SYS: a constant fitted on the total scale (K₁ and K₂ under `"Roy 1993"`, `"Lueker 2000"`, `"Sulpis 2020"`, `"SB 2020"` and `"Papadimitriou 2018"`; K_B; K_H2S; K_NH3 under `"Clegg"`) is moved to the seawater scale with the 1 atm conversion factor, corrected, and returned to the total scale with the factor at pressure. Applying the correction directly on the total scale instead shifts pH by about +1e-3 at 5000 dbar and +2e-3 at 11 000 dbar. Under KGen, K_H2S and a `"Clegg"` K_NH3 are still corrected on the total scale, because Kgen gives no surface K_SO4 and K_F; this matters only with `H2ST` or `NH4T` at depth.

#### KGen's composition correction — `MyAMI_mode`

Only `"approximate"`*, the polynomial approximation to MyAMI, is available in Kgen.jl; the full MyAMI model (`"calculate"`) is in the Python Kgen only, and asking for it is an error.

---

#### Non-Configurable Standard Constants

The following constants use single default parameterizations across the package:

* **$\text{CO}_2$ Solubility ($K_0$):** Weiss (1974)
* **Boric Acid ($K_\text{B}$):** Dickson (1990)
* **Water Dissociation ($K_\text{W}$):** Millero (1995)
* **Phosphoric Acid ($K_\text{P1}, K_\text{P2}, K_\text{P3}$) & Silicic Acid ($K_\text{Si}$):** Yao & Millero (1995)
* **Solubility Products ($K_\text{spA}, K_\text{spC}$):** Mucci (1983)
* **Hydrogen Sulfide ($K_\text{H2S}$):** Millero et al. (1988) / Yao & Millero (1995)

### Choosing the constants (`K_settings`)

Every choice of parameterisation is one field of `K_settings`. Give it as a NamedTuple naming only the choices that differ from the defaults, anywhere a calculation is set up:

```julia
using CarbonateCalculator

my_settings = (K1_K2_method = "Roy 1993", KF_method = "Perez")

# A single calculation
carbon_system(TA = 2300.0, DIC = 2000.0, K_settings = my_settings)

# A solver: the settings are checked once, here, not per sample
solve = CarbonateSystem(:carbon; varying = (:TA, :DIC), K_settings = my_settings)

# Constants on their own
CarbonateCalculator.calculate_constants(temp_c = 20.0, sal = 35.0, K_settings = my_settings)
```

`KSettings(...)` takes the same fields and gives the settings as a value you can inspect and reuse; `KSettings()` is the defaults, `KSettings(s; KF_method = "Dickson")` copies `s` with a change, and it prints as the call that rebuilds it. `?KSettings` in the REPL lists every option.

| Field | What it chooses | Options | Default |
| :--- | :--- | :--- | :--- |
| `K1_K2_method` | Carbonic acid $K_1, K_2$ | see the K₁/K₂ table above | `"KGen"` |
| `KSO4_method` | Bisulphate $K_\text{SO4}$ | `"Dickson"`, `"Khoo"`, `"WM13"` | `"Dickson"` |
| `KF_method` | Hydrogen fluoride $K_\text{F}$ | `"Dickson"`, `"Perez"` | `"Dickson"` |
| `KNH3_method` | Ammonium $K_\text{NH3}$ | `"Millero"`, `"Clegg"` | `"Millero"` |
| `BT_method` | Total boron from salinity | `"Uppstrom"`, `"Lee"`, `"KSK18"`, `"Culkin"`, `"freshwater"` | `"Uppstrom"` |
| `Ca_method` | Calcium from salinity | `"modern"`, `"Culkin"`, `"RT67"` | `"modern"` |
| `fH_method` | NBS activity coefficient $f_\text{H}$ | `"Takahashi 1982"`, `"Peng 1987"`, `"freshwater"` | `"Takahashi 1982"` |
| `pressure_method` | Pressure-correction coefficients | `"Millero 1995"`, `"Millero 1983"` | `"Millero 1995"` |
| `MyAMI_mode` | KGen's Ca/Mg correction | `"approximate"` | `"approximate"` |

A supplied `BT`, `Ca` or `Mg` value always overrides `BT_method` or `Ca_method`.

Each field is independent: no choice changes another. The K₁/K₂ parameterisations that used to switch other constants silently (GEOSECS, Peng et al., freshwater) now warn instead — see [Companion sets](#companion-sets).

### Warnings

The package warns, and still computes, in these cases. Each warning is shown once per distinct cause, so a broadcast over thousands of samples gives one warning for each way the data falls outside a range rather than one per sample (and so cannot say how many samples were affected).

When the settings are made (building a solver, or a single calculation):

* **Missing companions.** A K₁/K₂ parameterisation with [companion sets](#companion-sets), selected without any of them: *"You've chosen "Mehrbach 1973" for K1_K2_method. This K1/K2 method is intended to be used with one of these sets of companion constants: … You currently have "Takahashi 1982" selected for fH_method, "Uppstrom" selected for BT_method, "modern" selected for Ca_method."*
* **Companion used alone.** A companion-only option with a different K₁/K₂ parameterisation, e.g. `fH_method = "Peng 1987"` with `"Lueker 2000"`.
* **Ignored under KGen.** `KSO4_method`, `KF_method` or `pressure_method` chosen alongside `K1_K2_method = "KGen"`.

Per calculation (including [`adjust`](#changing-temperature--pressure-adjust), which checks the new conditions):

* **Out of range.** The water's temperature or salinity is outside the range the K₁/K₂ parameterisation was fitted over: *"You currently have selected "Lueker 2000" for K1_K2_method, however this method is not recommended for temp_c = 40.0 °C: it was fitted over 2.0 to 35.0 °C (Lueker et al. (2000); CO2SYS option 10)."*
* **Non-modern composition.** A `Ca` or `Mg` other than modern seawater with any K₁/K₂ parameterisation but KGen, which is the only one that corrects the constants for seawater composition.

Neither check runs when ready-made constants are passed with `Ks`, since they were not computed from the settings. To silence a warning you have judged harmless, use Julia's logging, e.g. `with_logger(NullLogger()) do … end` from `Logging`.

### Moving from `K_method`

The flat keywords `K_method`, `KSO4_method`, `KF_method`, `KNH3_method`, `BT_method`, `Ca_method` and `MyAMI_mode` still work, with a warning that prints the `K_settings` to write instead. The old K₁/K₂ names that switched other constants are translated into those constants explicitly, so old scripts give the same numbers:

| Old | New `K_settings` |
| :--- | :--- |
| `K_method = "Lueker 2000"`, or any other name not listed here | `(K1_K2_method = "Lueker 2000",)` |
| `K_method = "Mehrbach 1973 A"` | `(K1_K2_method = "Mehrbach 1973", BT_method = "Culkin")` |
| `K_method = "Mehrbach 1973 B"` | `(K1_K2_method = "Mehrbach 1973", fH_method = "Peng 1987")` |
| `K_method = "Millero 1979"` | `(K1_K2_method = "Millero 1979", fH_method = "freshwater", BT_method = "freshwater", pressure_method = "Millero 1983")` |
| `K_method = "automatic"` | removed — choose a `K1_K2_method`; an unsuitable one now warns per calculation |
| `KF_method = "Perez"` etc. | `(KF_method = "Perez",)` etc. |

The old Mehrbach translations reproduce what the package did, which was a subset of each CO2SYS set: Mehrbach 1973 B did not use Culkin boron or calcium, for instance. The companion warning will point that out; the full sets are in [Companion sets](#companion-sets). Under `"Millero 1979"`, K_B and the phosphate and silicate constants are no longer zeroed, which only mattered with a non-zero phosphate total, where it gave the wrong sign for phosphate alkalinity.

Inside `K_settings`, `"Mehrbach 1973 A"`, `"Mehrbach 1973 B"` and `"automatic"` are errors that name the replacement, and so is `K_method` as a field (it is now `K1_K2_method`). Results record `res.settings.K_settings` instead of the flat names.

## References

All equilibrium constant parameterizations, chemical formulations, and isotopic models implemented in `CarbonateCalculator.jl` are compiled in [`references.bib`](references.bib). 

If you use specific constant selections in your work (e.g., Lueker et al., 2000; Lee et al., 2010), you can import `references.bib` directly into Zotero, Mendeley, or BibTeX to cite the original oceanographic literature.