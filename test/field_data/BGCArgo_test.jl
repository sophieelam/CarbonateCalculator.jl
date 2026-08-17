# SOCCOM BGC-Argo floats — autonomous profiles from the Southern Ocean, served by NOAA
# PolarWatch's ERDDAP.
#
# pH is measured in situ by an ISFET sensor; TA and DIC are the LIAR/LIN empirical estimates
# distributed with the float data. So this checks the package against a measured pH at
# genuine in-situ temperature *and* pressure — the only dataset here that does — but the
# tolerance has to absorb the LIAR estimate's own error (~5-10 µmol/kg on TA), which is
# larger than anything the calculator contributes.

using Test
using CSV, DataFrames, Statistics, Printf, Downloads
using CarbonateCalculator

const BGCARGO_FILE = "SOCCOM_BGC_Argo.csv"

# One week of profiles. ERDDAP serves this as CSV with a units row beneath the header.
#
# ⚠️ This URL is dead. As of the last check the PolarWatch ERDDAP returns 404 for
# `SOCCOM_BGC_Argo`, and searching that server for "SOCCOM" or "BGC Argo" returns nothing —
# the dataset has been retired, not merely re-indexed, so adjusting the time range will not
# revive it. `with_dataset` therefore skips this file with a warning rather than failing.
# Replacing it means picking a new source (the SOCCOM data portal, or an Argo GDAC
# synthetic profile index) and remapping the columns; the checks below need only in-situ
# temperature, pressure, salinity, measured pH, and TA/DIC estimates.
const BGCARGO_URL = "https://polarwatch.noaa.gov/erddap/tabledap/SOCCOM_BGC_Argo.csv" *
    "?depth%2Cpressure%2Ctemperature%2Csalinity%2CpH_insitu%2CTALK_LIAR%2CDIC_LIAR" *
    "&time%3E=2021-10-15T00%3A00%3A00Z&time%3C=2021-10-22T14%3A47%3A00Z"

# `Downloads` is a stdlib, so this needs no HTTP.jl dependency the way the original did.
fetch_BGCArgo(path) = Downloads.download(BGCARGO_URL, path)


"Load the ERDDAP extract and reduce it to physically plausible rows."
function load_BGCArgo(path)
    raw = CSV.read(path, DataFrame, header = 1, skipto = 3)   # row 2 is units
    nrow(raw) == 0 && return raw

    df = DataFrame(
        depth = Float64.(raw[:, 1]),
        pres_bar = Float64.(raw[:, 2]),
        temp_c = Float64.(raw[:, 3]),
        sal = Float64.(raw[:, 4]),
        pH = Float64.(raw[:, 5]),
        TA = Float64.(raw[:, 6]),
        DIC = Float64.(raw[:, 7]),
    )

    # Floats report -999 for bad data, which would otherwise sail through as a number.
    df.temp_c = [(-5 < x < 40) ? x : missing for x in df.temp_c]
    df.sal = [(0 < x < 50) ? x : missing for x in df.sal]
    df.pres_bar = [(0 <= x < 6000) ? x : missing for x in df.pres_bar]
    df.pH = [(6 < x < 9) ? x : missing for x in df.pH]
    df.TA = [(1000 < x < 3000) ? x : missing for x in df.TA]
    df.DIC = [(1000 < x < 3000) ? x : missing for x in df.DIC]

    dropmissing!(df, [:pres_bar, :temp_c, :sal, :pH, :TA, :DIC])
    df.pres_bar .= df.pres_bar ./ 10.0   # dbar to bar

    return df
end


"""
Compute pH from the LIAR TA and DIC estimates, at the float's in-situ conditions.

Only pH is checked against a measurement. Computing TA from pH+DIC (and vice versa) would
be comparing the package against LIAR's regression rather than against anything measured,
which is why the original script's three-way comparison is reduced to one here.
"""
function solve_BGCArgo(df)
    pH = fill(NaN, nrow(df))

    for i in 1:nrow(df)
        try
            # No `Ks=`: see the note in GLODAP_test.jl.
            pH[i] = carbon_system(TA = df.TA[i], DIC = df.DIC[i], temp_c = df.temp_c[i],
                                  sal = df.sal[i], pres_bar = df.pres_bar[i],
                                  unit = "umol").pHtot
        catch
        end
    end

    return pH
end


@testset "SOCCOM BGC-Argo" begin
    with_dataset("SOCCOM BGC-Argo", BGCARGO_FILE, fetch_BGCArgo) do path
        df = load_BGCArgo(path)

        if nrow(df) == 0
            @warn "BGC-Argo extract has no usable rows — skipping. See the note on " *
                  "BGCARGO_URL: the SOCCOM dataset has been retired from PolarWatch and " *
                  "needs replacing with a live source."
            return
        end

        @info "BGC-Argo: $(nrow(df)) profiles points, $(round(minimum(df.temp_c), digits=1)) " *
              "to $(round(maximum(df.temp_c), digits=1)) °C, to " *
              "$(round(maximum(df.pres_bar) * 10)) dbar"

        computed = solve_BGCArgo(df)

        # Tolerance is set by the LIAR TA/DIC estimates, not by the calculator: a ~6 µmol/kg
        # TA error alone moves computed pH by ~0.01.
        check_agreement("pH from LIAR TA and DIC", df.pH, computed;
                        max_median = 0.05, max_iqr = 0.05)

        if FIELD_PLOTS
            figures = joinpath(FIGURE_DIR, "BGCArgo")
            mkpath(figures)
            savefig(comparison_figure(df.pH, computed, "pH", "Depth", df.depth;
                                      dataset = "BGC-Argo"),
                    joinpath(figures, "BGCArgo_pH_Comparison.png"))
            @info "BGC-Argo: figures written to $figures"
        end
    end
end
