"""
    CarbonateCalculator

Marine carbonate, boron and boron-isotope chemistry.

Every calculation goes through one solver, [`CarbonateSystem`](@ref), whose *scope* names the
subsystems to compute. The functions below are presets over it:

| function | scope | solves |
|---|---|---|
| [`carbon_system`](@ref) | `(:carbon,)` | the carbonate system |
| [`whole_system`](@ref) | `(:carbon, :boron, :isotopes)` | all three together |
| [`boron_system`](@ref) | `(:boron, :isotopes)` | boron speciation and its isotopes |
| [`boron_isotopes`](@ref) | `(:isotopes,)` | the isotopes alone |

```julia
carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0).pHtot
```

For many samples, name the varying parameters once and broadcast — Julia cannot broadcast
over keyword arguments, so naming them up front is what makes a keyword API vectorisable:

```julia
solve = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c), K_method = "Lueker 2000")
results = solve.([2300.0, 2310.0, 2290.0], [2000.0, 2010.0, 1990.0], [25.0, 10.0, 2.0])
```

The vectors are typically columns of a table. Anything varying per sample goes in `varying`;
anything constant stays a keyword, stated once and unable to drift between rows. Names are
checked when the solver is built, so a mistake surfaces there rather than on row 700,000.

# Uncertainties

Name the uncertain parameters in `varying_errors` and give a σ for each, positionally after the
values. `result.err` then carries a matching uncertainty for every computed quantity:

```julia
solve = CarbonateSystem(:carbon; varying = (:TA, :DIC), varying_errors = (:TA, :DIC),
                        temp_c = 20.0)
solve(2300.0, 2000.0, 2.0, 2.0).err.pHtot
```

Uncertainties are propagated by automatic differentiation, so every derived quantity gets one
without anyone writing an error formula. A σ shared by every sample is just a number, since
scalars broadcast against vectors:

```julia
solve.([2300.0, 2310.0, 2290.0], [2000.0, 2010.0, 1990.0], 2.0, 2.0)
```

Inputs are assumed uncorrelated. Where one input reaches the answer by several routes — salinity
enters every equilibrium constant — the correlation is handled exactly, because it is one input
differentiated once.

# Collection conditions

A sample is usually measured on deck and reported at the depth it came from.
[`at_collection_conditions`](@ref) re-solves a result at another temperature and pressure,
carrying the conservative totals, the seawater composition and every method choice across, so
nothing has to be restated:

```julia
measured = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0)
at_collection_conditions(measured; temp_c = 2.0, pres_bar = 400.0).pHtot
```

It broadcasts too, taking the same arguments positionally, and `σ_temp_c` / `σ_pres_bar` add
uncertainty in the collection conditions — which are often less well known than the conditions
the sample was measured at:

```julia
at_collection_conditions.([measured, measured], [2.0, 4.0], 400.0, [0.5, 0.2])
```

# Results in a table

A vector of results is a Tables.jl table, one row per sample, so it goes straight into a
`DataFrame` (or any other Tables.jl sink):

```julia
using DataFrames
DataFrame(solve.([2300.0, 2310.0], [2000.0, 1990.0], 2.0, 2.0))
```

Columns are the computed values, and where a result carries uncertainties each value gains a
`σ_` column beside it — `pHtot` and `σ_pHtot`. Results that do not all carry the same columns
are filled with `nothing` rather than silently truncated.

The source is laid out by role: `parameters.jl` declares what the inputs are, `pipeline.jl`
how a system is solved, `system.jl` the solver itself, and `result.jl`, `validation.jl`,
`conditions.jl` and `display.jl` the result type, the input checks, the two-stage conditions
API and printing.
"""
module CarbonateCalculator

# Chemistry: the speciation equations and the equilibrium constants.
include("carbon.jl")  # all machinery for aqueous carbonate speciation and the carbonate system
using .Carbon
include("boron.jl")  # all machinery for aqueous boron speciation
using .Boron
include("boron_isotopes.jl")  # all machinery for boron isotopes and their fractionation
using .Isotopes
include("constants.jl")  # all constants used in calculations, including the equilibrium constants and their temperature/salinity/pressure dependencies, as well as salinity-derived seawater composition.
using .Constants
include("helpers.jl")  # currently just a helper function to convert pH scales (may be broken?)
using .Helpers

using Printf
using ForwardDiff, LinearAlgebra, Roots

# The package proper, in dependency order.
include("errors.jl")        # Machinery for uncertainty propagation
include("result.jl")        # CarbonateResult data type, which carries the results of a calculation and its uncertainties
include("tables.jl")        # Tables.jl interface to a CarbonateResult (used to convert a vector of results to a DataFrame)
include("parameters.jl")    # A map of possible inputs and degrees of freedom, defining calculation scopes for the solvers
include("pipeline.jl")      # All the machinery that the solver uses to do the calculations
include("validation.jl")    # Checks for validity of inputs
include("display.jl")       # A pretty-printing method for CarbonateResult, and the logic to decide which variables to print
include("system.jl")        # CarbonateSystem, the function to generate solvers, and some user-facing presets.
include("conditions.jl")    # at_collection_conditions function, used to re-solve a result at the conditions the sample was collected at, rather than the conditions it was measured at.
include("gradients.jl")     # Functions for calculating buffer factors and gradients (e.g. δpCO₂/δTA, δpH/δDIC, etc.) and their uncertainties

export CarbonateSystem,  # the main solver
       carbon_system, whole_system, boron_system, boron_isotopes,  # user-facing convenience functions
       propagate_errors, 
       at_collection_conditions,  # re-solve a result at the conditions the sample was collected at
       calc_gradient, calc_relative_gradient, calc_gradient_uncertainty,  # buffer factors
       revelle_factor, with_gradient

end # module
