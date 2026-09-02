using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Revise
using CarbonateCalculator

using DataFrames

## Basic usage example

res = carbon_system(TA=2300, DIC=2100)

df = DataFrame([res])

df.pHtot

## Advanced usage example

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

## Case study: a palaeo use case, reconstructing pCO2 from δ¹¹B of borate, assuming constant seawater composition and isotopic composition.

palaeo_solver = CarbonateSystem(
    :carbon, :isotopes,  # we need both the carbon and boron isotope machinery
    varying=(:δBOH₄, :temp_c),
    varying_errors=(:δBOH₄,),
    sal=35.0, TA=2300  # constant seawater composition
)

# get a range of 50 d11B values for borate, and their errors, to reconstruct pCO2 from
d11B_values = range(15, 18, length=50)
d11B_errors = 0.5 .+ 0.1 .* randn(50)
temp_c = 20.0 .+ 2.0 .* randn(50)  # add some random temperature variation to the palaeo reconstruction

# run the solver over the range of values, and get the results
palaeo_results = palaeo_solver.(d11B_values, temp_c, d11B_errors)

df = DataFrame(palaeo_results)

using Plots

plot(df.pHtot,ribbon=df.σ_pHtot, label="pH ± σ", xlabel="Sample index", ylabel="pH")
