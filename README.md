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
| `BT` | Total boron | `unit` | from salinity (`BT_method`) | all |
| `ST` | Total sulphate | mol/kg | from salinity | all |
| `FT` | Total fluoride | mol/kg | from salinity | all |
| `PT`, `SiT` | Total phosphate, silicate | `unit` | 0.0 | carbon |
| `H2ST`, `NH4T` | Total hydrogen sulphide, ammonia | `unit` | 0.0 | carbon |
| `Ca`, `Mg` | Calcium, magnesium, given at S = 35 | mol/kg | modern seawater (`Ca_method`) | all |
| **Settings** | | | | |
| `unit` | Concentration unit: `"mol"`, `"mmol"`, `"umol"`, `"nmol"`, `"pmol"`, `"fmol"` (per kg) | — | `"umol"` | all |
| `K_method` | K₁/K₂ parameterisation; `"KGen"` sets every constant | — | `"KGen"` | all |
| `KSO4_method`, `KF_method`, `KNH3_method` | Bisulphate, HF and ammonium constants (no effect with KGen, except `KNH3_method`) | — | see [Chemical Constants](#chemical-constants) | all |
| `BT_method`, `Ca_method` | Total boron and calcium from salinity | — | `"Uppstrom"`, `"modern"` | all |
| `MyAMI_mode` | KGen's Ca/Mg correction: `"approximate"` or `"calculate"` | — | `"approximate"` | all |
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

Besides these values, a result also carries `res.Ks` (the equilibrium constants used, in mol/kg on the total pH scale, with `KS` and `KF` on the free scale), `res.unit`, and `res.settings` (every setting above). `revelle_factor(res)` and `calc_gradient` give buffer factors.

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

Note that by default `K_method="Kgen"`, where *ALL* K values are calculated using [Kgen.jl](https://palaeocarbonatechemistry.github.io/Kgen/), overriding K-specific method flags below.
[**KGen**](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GC011417) provides 'best practice' K values for modern seawater, and provides a polynomial parameterisation for the influence of Mg and Ca variation in palaeo-seawater. 

Equilibrium constant parameterizations and seawater composition models can be specified using `*_method` arguments when constructing a solver or calling preset functions. Full BibTeX entries for all supported methods are available in [`references.bib`](references.bib).

#### Carbonate ($K_1, K_2$) — `K_method`

| Key / Method | Reference | Temp (°C) | Salinity | pH Scale & Notes |
| :--- | :--- | :--- | :--- | :--- |
| `"KGen"`* | Whiteford et al. (2025) | -2 to 50 | 0 to 50 | Total scale; modern & palaeo seawater |
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


#### Sulfate ($K_\text{SO4}$, $K_\text{HSO4}$) — `KSO4_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Dickson"`* | Dickson (1990) | Standard oceanographic formulation. |
| `"Kuo"` | Khoo et al. (1977) | Free pH scale formulation. |
| `"WM13"` | Waters & Millero (2013) | Extended temperature and salinity parameterization.

**This setting has no effect if `K_method="Kgen"`.**

#### Hydrogen Fluoride ($K_\text{F}$) — `KF_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Dickson"`* | Dickson & Riley (1979) | Standard oceanographic formulation. |
| `"Perez"` | Pérez & Fraga (1987) | Alternative fit ($S = 10\text{--}40$, $T = 9\text{--}33^\circ\text{C}$).

**This setting has no effect if `K_method="Kgen"`.**


#### Ammonia ($K_\text{NH3}$) — `KNH3_method`

| Key / Method | Reference | Formulation / Description |
| :--- | :--- | :--- |
| `"Millero"`* | Yao & Millero (1995) | Seawater pH scale formulation. |
| `"Clegg"` | Clegg & Whitfield (1995) | Total pH scale fit ($S = 0\text{--}40$, $T = -2\text{--}40^\circ\text{C}$).

---

#### Non-Configurable Standard Constants

The following constants use single default parameterizations across the package:

* **$\text{CO}_2$ Solubility ($K_0$):** Weiss (1974)
* **Boric Acid ($K_\text{B}$):** Dickson (1990)
* **Water Dissociation ($K_\text{W}$):** Millero (1995)
* **Phosphoric Acid ($K_\text{P1}, K_\text{P2}, K_\text{P3}$) & Silicic Acid ($K_\text{Si}$):** Yao & Millero (1995)
* **Solubility Products ($K_\text{spA}, K_\text{spC}$):** Mucci (1983)
* **Hydrogen Sulfide ($K_\text{H2S}$):** Millero et al. (1988) / Yao & Millero (1995)

## References

All equilibrium constant parameterizations, chemical formulations, and isotopic models implemented in `CarbonateCalculator.jl` are compiled in [`references.bib`](references.bib). 

If you use specific constant selections in your work (e.g., Lueker et al., 2000; Lee et al., 2010), you can import `references.bib` directly into Zotero, Mendeley, or BibTeX to cite the original oceanographic literature.