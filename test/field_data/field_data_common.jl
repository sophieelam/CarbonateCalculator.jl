# Shared scaffolding for block 3 of the suite: accuracy against real-world measurements.
#
# Each dataset in this directory has the carbonate system over-determined — it reports more
# parameters than are needed to solve it — so a computed parameter can be compared against
# the measured one. That is the only thing in the suite that says anything about accuracy
# rather than self-consistency.
#
# These files used to be scripts that downloaded data, ran the calculator and wrote PNGs
# without a single assertion. Every row's solve sat in a bare `catch e → NaN`, so when the
# `Ks=` argument was wrong the calculator raised on *every* row and the only symptom was an
# empty plot nobody was looking at. Hence `min_solved` below: the assertion that the
# calculation ran at all matters as much as the assertion that it was accurate.

module FieldData

using Test, Statistics, Printf

export MAKE_PLOTS, DATA_DIR, FIGURE_DIR, dataset_path, with_dataset, agreement,
       check_agreement, optional_using

"""
    optional_using(name) -> Bool

Load a package that only block 3 needs, reporting how to install it if it is absent.

Kept out of `test/Project.toml` deliberately. Blocks 1 and 2 need nothing beyond `CSV`,
`DataFrames` and the stdlibs, and making everyone who runs `Pkg.test()` install `Plots` to
satisfy a block that is off by default is the wrong trade. Block 3 degrades to a warning
instead.
"""
function optional_using(name::Symbol)
    try
        @eval Main using $name
        return true
    catch
        @warn "Optional dependency $name is not installed, so part of block 3 will be " *
              "skipped. Install it with:\n" *
              "    julia --project=test -e 'using Pkg; Pkg.add(\"$name\")'"
        return false
    end
end

"""
Whether to write comparison figures as well as checking the numbers. Off by default:
plotting pulls in `Plots`, which is slow to load, and the figures are a benchmarking aid
rather than a check. Enable with `CC_PLOTS=true`.
"""
const MAKE_PLOTS = get(ENV, "CC_PLOTS", "false") == "true"

"""
Where downloaded datasets are cached. Persistent and gitignored (bar the committed GLODAP
subset), so the several hundred MB block 3 needs is fetched once per machine and every
later run reads from disk.
"""
const DATA_DIR = joinpath(@__DIR__, "data")

const FIGURE_DIR = joinpath(@__DIR__, "figures")

"Path a dataset is cached at. Everything in block 3 goes through here, so there is one
place that decides where downloads live."
dataset_path(filename::AbstractString) = joinpath(DATA_DIR, filename)

"""
    with_dataset(body, name, filename, fetch)

Run `body(path)` with the dataset cached at `dataset_path(filename)`.

`fetch(path)` is called **only when the file is not already on disk**, so the download
happens once per machine and re-runs are offline. If it fails, warn and return `nothing`
without running `body`: a dataset being unreachable means the network or the host is down,
not that the package has regressed, and block 3 should not turn a NOAA outage into a red
suite. The opposite convention is what makes people stop running the block at all.
"""
function with_dataset(body, name::AbstractString, filename::AbstractString, fetch)
    path = dataset_path(filename)

    if isfile(path)
        @info "$name: using cached $(basename(path)) ($(round(filesize(path) / 1e6, digits=1)) MB)"
    else
        mkpath(DATA_DIR)
        try
            @info "$name: not cached — downloading once to $path"
            fetch(path)
        catch err
            @warn "Could not fetch $name — skipping. This is a data availability " *
                  "problem, not a package regression." exception = err
            return nothing
        end
        isfile(path) || (@warn "$name still missing after download — skipping."; return nothing)
    end

    return body(path)
end

"""
    agreement(measured, computed) -> NamedTuple

Robust summary of `computed .- measured`.

Median and IQR rather than mean and standard deviation: field data carries occasional
badly-flagged rows, and a single outlier moves a mean enough to make any tolerance either
meaningless or unpassable. `solved_fraction` counts rows that produced a number at all.
"""
function agreement(measured, computed)
    difference = computed .- measured
    solved = .!isnan.(difference)
    valid = difference[solved]

    isempty(valid) && return (n = length(difference), n_solved = 0, solved_fraction = 0.0,
                              median_offset = NaN, iqr = NaN, p95 = NaN, worst = NaN)

    return (
        n = length(difference),
        n_solved = length(valid),
        solved_fraction = length(valid) / length(difference),
        median_offset = median(valid),
        iqr = quantile(valid, 0.75) - quantile(valid, 0.25),
        p95 = quantile(abs.(valid), 0.95),
        worst = maximum(abs.(valid)),
    )
end

"""
    check_agreement(label, measured, computed; max_median, max_iqr, min_solved=0.95)

Assert that the package reproduces a measured field parameter, and report the statistics.

Tolerances are set from observed agreement with headroom, not from what would be nice — see
`CLEANUP.md` for the measured values behind each. `min_solved` defaults to 0.95 rather than
1.0 because real data contains rows the solver legitimately cannot converge on.
"""
function check_agreement(label, measured, computed; max_median, max_iqr, min_solved = 0.95)
    stats = agreement(measured, computed)

    @printf("  %-24s solved %5d/%-5d (%5.1f%%)  median %+9.4g  IQR %8.4g  p95 %8.4g\n",
            label, stats.n_solved, stats.n, 100 * stats.solved_fraction,
            stats.median_offset, stats.iqr, stats.p95)

    @testset "$label" begin
        # The calculation ran at all. This is the check the `Ks=` bug would have failed.
        @test stats.solved_fraction >= min_solved
        # ...and was accurate.
        @test abs(stats.median_offset) <= max_median
        @test stats.iqr <= max_iqr
    end

    return stats
end

end # module
