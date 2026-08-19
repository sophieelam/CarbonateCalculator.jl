# CODAP-NA v2021 — North American coastal ocean data, discrete bottle samples.
#
# This is the one dataset that validates the **two-stage conditions API** against real
# measurements, because CODAP reports the same sample in both states:
#
#   * `pH_TS_measured` at `TEMP_pH` and atmospheric pressure — the bench measurement;
#   * `pH_TS_insitu_measured` at `CTDTEMP_ITS90` and `CTDPRES` — the same water where it
#     was collected.
#
# So `carbon_system` is checked at the conditions the sample was measured at, and
# `at_collection_conditions` is checked at the conditions it was collected at,
# against a published value for each. Nothing in block 1 can test that: it can only show
# the two stages agree with each other.
#
# Coastal rather than open-ocean, so it also exercises a wider salinity range and higher
# nutrient concentrations than GLODAPv2, where the phosphate and silicate alkalinity
# contributions stop being negligible.

using Test
using CSV, DataFrames, Statistics, Printf, Downloads
using CarbonateCalculator

const CODAP_FILE = "CODAP_NA_v2021.csv"
const CODAP_URL = "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0219960/" *
                  "Data_CODAP/CSV/CODAP_NA_v2021.csv"

fetch_CODAP(path) = Downloads.download(CODAP_URL, path)

# Tolerances, set from measured agreement with headroom — see CLEANUP.md. The in-situ
# comparison is the looser of the two because CODAP's own in-situ pH is itself a
# temperature/pressure extrapolation of the bench measurement, so the two sides disagree
# about the extrapolation as well as about the chemistry.
const CODAP_LAB_MAX_MEDIAN = 0.02
const CODAP_LAB_MAX_IQR = 0.04
const CODAP_INSITU_MAX_MEDIAN = 0.02
const CODAP_INSITU_MAX_IQR = 0.05


"Load CODAP-NA and reduce it to rows carrying both the lab and the in-situ state."
function load_CODAP(path)
    df = CSV.read(path, DataFrame, missingstring = ["-999", "-9999", "NaN"])

    for (from, to) in ("TEMP_pH" => :lab_temp_c, "CTDTEMP_ITS90" => :insitu_temp_c,
                       "CTDPRES" => :insitu_pres_bar,
                       "recommended_Salinity_PSS78" => :sal,
                       "DIC" => :DIC, "TALK" => :TA,
                       "pH_TS_measured" => :lab_pH,
                       "pH_TS_insitu_measured" => :insitu_pH,
                       "Silicate" => :SiT, "Phosphate" => :PT)
        hasproperty(df, Symbol(from)) && rename!(df, Symbol(from) => to)
    end

    # CODAP mixes numeric and string columns depending on the export, so parse defensively.
    parse_number(x) = ismissing(x) ? missing : something(tryparse(Float64, string(x)), missing)
    for column in (:lab_temp_c, :insitu_temp_c, :insitu_pres_bar, :sal, :DIC, :TA,
                   :lab_pH, :insitu_pH, :SiT, :PT)
        hasproperty(df, column) && (df[!, column] = parse_number.(df[!, column]))
    end

    # A missing nutrient means "not measured", which for alkalinity purposes is zero.
    for nutrient in (:SiT, :PT)
        df[!, nutrient] = hasproperty(df, nutrient) ? coalesce.(df[!, nutrient], 0.0) :
                                                      zeros(nrow(df))
    end

    dropmissing!(df, [:DIC, :TA, :lab_pH, :insitu_pH, :lab_temp_c, :insitu_temp_c,
                      :insitu_pres_bar, :sal])
    filter!(row -> row.sal > 20 && row.DIC > 1000 && row.TA > 1000, df)
    df.insitu_pres_bar .= df.insitu_pres_bar ./ 10.0   # dbar to bar

    return df
end


"""
Solve each sample at the bench, then carry it to the conditions it was collected at.

Stage 2 goes through `at_collection_conditions` rather than solving directly at the
in-situ conditions. The two are equivalent — TA and DIC are conservative — which is exactly
why this is worth doing against real data: it checks the machinery that carries the totals,
the seawater composition and every method choice across the condition change, and compares
the answer to a published in-situ pH rather than to another calculation of our own.
"""
function solve_CODAP(df)
    lab_pH = fill(NaN, nrow(df))
    insitu_pH = fill(NaN, nrow(df))

    for i in 1:nrow(df)
        try
            # No `Ks=`: see the note in GLODAP_test.jl. `pres_bar = 0` because this is a
            # bench measurement on a collected sample — the sample is no longer at depth.
            measured = carbon_system(TA = df.TA[i], DIC = df.DIC[i],
                                     temp_c = df.lab_temp_c[i], sal = df.sal[i],
                                     pres_bar = 0.0, PT = df.PT[i], SiT = df.SiT[i],
                                     unit = "umol")
            lab_pH[i] = measured.pHtot

            collected = at_collection_conditions(measured;
                                                         temp_c = df.insitu_temp_c[i],
                                                         pres_bar = df.insitu_pres_bar[i])
            insitu_pH[i] = collected.pHtot
        catch
        end
    end

    return (lab = lab_pH, insitu = insitu_pH)
end


@testset "CODAP-NA v2021" begin
    with_dataset("CODAP-NA v2021", CODAP_FILE, fetch_CODAP) do path
        df = load_CODAP(path)
        @info "CODAP-NA: $(nrow(df)) samples; bench " *
              "$(round(minimum(df.lab_temp_c), digits=1))-$(round(maximum(df.lab_temp_c), digits=1)) °C, " *
              "in situ $(round(minimum(df.insitu_temp_c), digits=1))-" *
              "$(round(maximum(df.insitu_temp_c), digits=1)) °C to " *
              "$(round(maximum(df.insitu_pres_bar) * 10)) dbar"

        computed = solve_CODAP(df)

        check_agreement("pH at bench conditions", df.lab_pH, computed.lab;
                        max_median = CODAP_LAB_MAX_MEDIAN, max_iqr = CODAP_LAB_MAX_IQR)
        check_agreement("pH recalculated in situ", df.insitu_pH, computed.insitu;
                        max_median = CODAP_INSITU_MAX_MEDIAN,
                        max_iqr = CODAP_INSITU_MAX_IQR)

        if FIELD_PLOTS
            figures = joinpath(FIGURE_DIR, "CODAP")
            mkpath(figures)
            savefig(comparison_figure(df.lab_pH, computed.lab, "pH (bench)", "Silicate",
                                      df.SiT; dataset = "CODAP-NA"),
                    joinpath(figures, "CODAP_pH_bench.png"))
            savefig(comparison_figure(df.insitu_pH, computed.insitu, "pH (in situ)",
                                      "Pressure", df.insitu_pres_bar; dataset = "CODAP-NA"),
                    joinpath(figures, "CODAP_pH_insitu.png"))
            @info "CODAP-NA: figures written to $figures"
        end
    end
end
