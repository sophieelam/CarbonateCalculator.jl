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
solve = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c, :sal),
                        K_method = "Lueker 2000")
results = solve.(df.talk, df.tco2, df.temperature, df.salinity)
```

Anything varying per sample goes in `varying`; anything constant stays a keyword, stated once
and unable to drift between rows. Names are checked when the solver is built, so a mistake
surfaces there rather than on row 700,000.

To report a sample at the conditions it was *collected* at rather than measured at, pass a
result to [`at_collection_conditions`](@ref).

The source is laid out by role: `parameters.jl` declares what the inputs are, `pipeline.jl`
how a system is solved, `system.jl` the solver itself, and `result.jl`, `validation.jl`,
`conditions.jl` and `display.jl` the result type, the input checks, the two-stage conditions
API and printing.
"""
module CarbonateCalculator

# Chemistry: the speciation equations and the equilibrium constants.
include("carbon.jl")
using .Carbon
include("boron.jl")
using .Boron
include("boron_isotopes.jl")
using .Isotopes
include("constants.jl")
using .Constants
include("helpers.jl")
using .Helpers

using Printf
using ForwardDiff, LinearAlgebra, Roots

# The package proper, in dependency order.
include("errors.jl")        # uncertainty propagation
include("result.jl")        # CarbonateResult
include("tables.jl")        # a vector of results is a Tables.jl table
include("parameters.jl")    # what the inputs are
include("pipeline.jl")      # how a system is solved
include("validation.jl")    # input checks
include("display.jl")       # how a result prints
include("system.jl")        # the solver, and the presets over it
include("conditions.jl")    # re-solving at other conditions

export CarbonateSystem,  # the main solver
       carbon_system, whole_system, boron_system, boron_isotopes,  # user-facing convenience functions
       propagate_errors, 
       at_collection_conditions  # re-solve a result at the conditions the sample was collected at

end # module
