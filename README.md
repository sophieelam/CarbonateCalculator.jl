# CarbonateCalculator.jl

[![Run Tests and Coverage](https://github.com/sophieelam/CarbonateCalculator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sophieelam/CarbonateCalculator.jl/actions/workflows/ci.yml)

[![codecov](https://codecov.io/github/sophieelam/carbonatecalculator.jl/graph/badge.svg?token=NH92EBMYRY)](https://codecov.io/github/sophieelam/carbonatecalculator.jl)

Welcome to CarbonateCalculator.jl!

## Installation

Equilibrium constants for palaeo seawater come from [Kgen](https://github.com/oscarbranson/Kgen), which is not
yet in the General registry, so it has to be added by URL first:

```julia
using Pkg
Pkg.add(url="https://github.com/oscarbranson/Kgen", rev="julia", subdir="julia/Kgen.jl")
Pkg.develop(path="path/to/CarbonateCalculator.jl")
```

The committed `Manifest.toml` already pins Kgen to that branch, so `Pkg.instantiate()` in a
checkout of this repository is enough on its own. Once Kgen is registered, both steps
collapse into a plain `Pkg.add`.

No Python is required. Constants are calculated by the pure-Julia Kgen, which means
`ForwardDiff` differentiates through them and error propagation works on every method,
including MyAMI. Python is needed only to regenerate the PyCO2SYS comparison baseline in
`test/CO2SYS_tests/`, and is declared in `test/Project.toml` and `test/CondaPkg.toml`
rather than as a dependency of the package.






