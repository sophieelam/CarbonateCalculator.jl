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

```jldoctest overview
julia> carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0).pHtot
8.1217859325142
```

For many samples, name the varying parameters once and broadcast — Julia cannot broadcast
over keyword arguments, so naming them up front is what makes a keyword API vectorisable:

```jldoctest overview
julia> solve = CarbonateSystem(:carbon; varying = (:TA, :DIC, :temp_c),
                               K_method = "Lueker 2000");

julia> results = solve.([2300.0, 2310.0, 2290.0], [2000.0, 2010.0, 1990.0], [25.0, 10.0, 2.0]);

julia> getproperty.(results, :pHtot)
3-element Vector{Float64}:
 8.045886181605269
 8.276233808810545
 8.40949674883427
```

The vectors are typically columns of a table. Anything varying per sample goes in `varying`;
anything constant stays a keyword, stated once and unable to drift between rows. Names are
checked when the solver is built, so a mistake surfaces there rather than on row 700,000.

Uncertainties are propagated by automatic differentiation — name them in `varying_errors` and
`result.err` carries a σ for every computed quantity. [`at_collection_conditions`](@ref)
re-solves a result at the depth a sample came from, and a vector of results is a Tables.jl
table that goes straight into a `DataFrame`. All three are covered below.

# Extended help

# Uncertainties

Name the uncertain parameters in `varying_errors` and give a σ for each, positionally after the
values. `result.err` then carries a matching uncertainty for every computed quantity:

```jldoctest overview
julia> solve = CarbonateSystem(:carbon; varying = (:TA, :DIC), varying_errors = (:TA, :DIC),
                               temp_c = 20.0);

julia> solve(2300.0, 2000.0, 2.0, 2.0).err.pHtot
0.004620992184145955
```

Uncertainties are propagated by automatic differentiation, so every derived quantity gets one
without anyone writing an error formula. A σ shared by every sample is just a number, since
scalars broadcast against vectors:

```jldoctest overview
julia> shared_sigma = solve.([2300.0, 2310.0, 2290.0], [2000.0, 2010.0, 1990.0], 2.0, 2.0);

julia> getproperty.(getproperty.(shared_sigma, :err), :pHtot)
3-element Vector{Float64}:
 0.004620992184145955
 0.004614381664736751
 0.0046276345876582175
```

Inputs are assumed uncorrelated. Where one input reaches the answer by several routes — salinity
enters every equilibrium constant — the correlation is handled exactly, because it is one input
differentiated once.

# Collection conditions

A sample is usually measured on deck and reported at the depth it came from.
[`at_collection_conditions`](@ref) re-solves a result at another temperature and pressure,
carrying the conservative totals, the seawater composition and every method choice across, so
nothing has to be restated:

```jldoctest overview
julia> measured = carbon_system(TA = 2300.0, DIC = 2000.0, temp_c = 20.0);

julia> at_collection_conditions(measured; temp_c = 2.0, pres_bar = 400.0).pHtot
8.254156644905194
```

It broadcasts too, taking the same arguments positionally, and `σ_temp_c` / `σ_pres_bar` add
uncertainty in the collection conditions — which are often less well known than the conditions
the sample was measured at:

```jldoctest overview
julia> cast = at_collection_conditions.([measured, measured], [2.0, 4.0], 400.0, [0.5, 0.2]);

julia> getproperty.(cast, :pHtot)
2-element Vector{Float64}:
 8.254156644905194
 8.223138522847922
```

# Results in a table

A vector of results is a Tables.jl table, one row per sample, so it goes straight into a
`DataFrame` (or any other Tables.jl sink):

```jldoctest overview
julia> using DataFrames

julia> df = DataFrame(solve.([2300.0, 2310.0], [2000.0, 1990.0], 2.0, 2.0));

julia> df.pHtot
2-element Vector{Float64}:
 8.1217859325142
 8.15376174968104

julia> df.σ_pHtot
2-element Vector{Float64}:
 0.004620992184145955
 0.004437417288202379
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
include("carbon.jl")          # aqueous carbonate speciation and the carbonate system
using .Carbon
include("boron.jl")           # aqueous boron speciation
using .Boron
include("boron_isotopes.jl")  # boron isotopes and their fractionation
using .Isotopes
include("constants.jl")       # equilibrium constants, their T/S/P dependencies, and
using .Constants              # salinity-derived seawater composition
include("helpers.jl")         # pH-scale conversion
using .Helpers

using Printf
using ForwardDiff, LinearAlgebra, Roots

# The package proper, in dependency order.
include("errors.jl")        # uncertainty propagation by automatic differentiation
include("result.jl")        # the CarbonateResult type, carrying values and uncertainties
include("tables.jl")        # Tables.jl interface, so a vector of results is a table
include("parameters.jl")    # every accepted parameter, its default, and its scope
include("pipeline.jl")      # the stages of a solve
include("validation.jl")    # input checks
include("display.jl")       # pretty-printing a CarbonateResult
include("system.jl")        # CarbonateSystem, the solver, and the four presets
include("conditions.jl")    # at_collection_conditions: re-solve at the collection conditions
include("gradients.jl")     # buffer factors and gradients, and their uncertainties

export CarbonateSystem,  # the solver
       carbon_system, whole_system, boron_system, boron_isotopes,  # presets over it
       propagate_errors,
       at_collection_conditions,  # re-solve at the conditions the sample was collected at
       calc_gradient, calc_relative_gradient, calc_gradient_uncertainty,  # buffer factors
       revelle_factor, with_gradient

# Part of the API, but reached by qualified name rather than brought into a user's namespace:
# `CarbonateResult` is what every calculation returns, and `calculate_constants` is the way to
# a constant bundle on its own.
public CarbonateResult, calculate_constants

end # module
