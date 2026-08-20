# Runs every `jldoctest` block in the package's docstrings, checking that each example still
# runs and still produces the output it claims.
#
# `Documenter.doctest` works straight off the module, so this needs no `docs/` directory and no
# site build, and Documenter stays a test-only dependency. If a docs site is added later, the
# same blocks are what it would check.

using Documenter, Logging

# Every doctest block is evaluated in a fresh module, so the `using` that makes `carbon_system`
# and friends visible has to be attached to the package rather than written into each example.
DocMeta.setdocmeta!(CarbonateCalculator, :DocTestSetup, :(using CarbonateCalculator);
                    recursive = true)

# The last digits of a float can differ between platforms and libm versions, and this package
# reaches its answers through exp, log and sqrt. Comparison therefore stops at eight decimal
# places: the filter deletes anything beyond that from the expected and actual output alike.
# Shorter values are compared in full.
const DOCTEST_FILTERS = [r"(?<=\d\.\d{8})\d+"]

# Documenter narrates each build step at Info level, which is six lines of noise per run. The
# threshold is Warn rather than silence so that a failing doctest still prints its diff.
with_logger(SimpleLogger(stderr, Logging.Warn)) do
    doctest(CarbonateCalculator; manual = false, doctestfilters = DOCTEST_FILTERS)
end
