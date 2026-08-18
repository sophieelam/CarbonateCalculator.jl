# GLODAPv2 — global bottle data, quality-controlled, with pH, DIC and TA all measured.
#
# The system is over-determined, so each of the three can be computed from the other two
# and compared against its own measurement. Agreement is expected at the level of the
# GLODAP internal-consistency literature (~0.01-0.03 pH, 5-15 µmol/kg), not at machine
# precision: the discrepancy is dominated by measurement uncertainty and by the equilibrium
# constants, not by this package.

using Test
using CSV, DataFrames, Statistics, Printf, Downloads
using CarbonateCalculator


# Rows to solve. The file holds ~143k usable rows; 2000 sampled evenly across it spans the
# full range of conditions (-2.0 to 30.4 °C, salinity 24.9 to 37.3, 0 to 650 bar) while
# keeping the block to a few seconds.
const GLODAP_SAMPLE_SIZE = 2000

const GLODAP_SUBSET = "GLODAPv2_pH_DIC_ALK_subset.csv"


"""
Download the GLODAPv2 master file, keep the quality-flagged rows and the columns this
comparison needs, and write the subset to `path`.

Only called when the subset is absent — the subset is committed, so in practice this runs
on no one's machine. It is kept because it is the only record of how the subset was made.
"""
function fetch_GLODAP(path)
    zip_path = joinpath(dirname(path), "GLODAPv2_Merged_Master_File.csv.zip")
    url = "https://www.glodap.info/glodap_files/v2.2020/GLODAPv2.2020_Merged_Master_File.csv.zip"

    optional_using(:ZipFile) || error("GLODAPv2 subset is absent and ZipFile is not " *
                                     "installed, so it cannot be rebuilt.")

    isfile(zip_path) || Downloads.download(url, zip_path)

    archive = Main.ZipFile.Reader(zip_path)
    try
        entry = only(filter(f -> f.name == "GLODAPv2.2020_Merged_Master_File.csv",
                            archive.files))
        gd = CSV.read(read(entry), DataFrame,
                      missingstring = ["-9999", "-9999.0", "-999.9"])

        # Flag 2 is "good" in the GLODAP convention; anything else is not usable here.
        for (value, flag) in ((:phtsinsitutp, :phtsinsitutpf), (:tco2, :tco2f),
                              (:talk, :talkf), (:salinity, :salinityf),
                              (:phosphate, :phosphatef), (:silicate, :silicatef))
            gd[!, value] = ifelse.(gd[!, flag] .== 2, gd[!, value], missing)
        end

        dropmissing!(gd, [:phtsinsitutp, :tco2, :talk, :temperature, :salinity,
                          :pressure, :silicate, :phosphate])

        CSV.write(path, select(gd, [
            :phts25p0, :phtsinsitutp, :tco2, :talk, :temperature, :salinity,
            :cruise, :station, :cast, :year, :month, :day, :hour,
            :latitude, :longitude, :bottomdepth, :maxsampdepth, :bottle,
            :pressure, :depth, :theta, :silicate, :phosphate]))
    finally
        close(archive)
    end
    return nothing
end


"Load the subset and reduce it to the rows this comparison can use."
function load_GLODAP(path)
    gd = CSV.read(path, DataFrame)

    dropmissing!(gd, [:phtsinsitutp, :temperature, :salinity, :tco2, :talk, :pressure,
                      :phosphate, :silicate])
    filter!(row -> row.salinity > 20.0 && row.tco2 > 1000.0 && row.talk > 1000.0, gd)
    filter!(row -> row.cruise != 270, gd)   # cruise 270 is known-bad in this file
    gd.pressure .= gd.pressure ./ 10.0      # dbar to bar

    # Sample evenly across the file rather than taking the first N rows, which would be a
    # single cruise and a single region.
    n = min(nrow(gd), GLODAP_SAMPLE_SIZE)
    return gd[round.(Int, range(1, nrow(gd), length = n)), :]
end


"""
Solve each row three ways: pH from TA+DIC, TA from pH+DIC, DIC from pH+TA.

A row that fails to converge becomes `NaN` rather than aborting the sweep — real data
contains rows the solver cannot handle. The count is asserted on, so failures cannot hide
the way they did when this file swallowed them into a plot.
"""
function solve_GLODAP(gd)
    pH  = fill(NaN, nrow(gd))
    TA  = fill(NaN, nrow(gd))
    DIC = fill(NaN, nrow(gd))

    for i in 1:nrow(gd)
        # No `Ks=` here. This used to precompute the constants with `calculate_constants`
        # and pass `row_K_results.Ks` — the inner 15-constant bundle where the outer
        # environment tuple is required — which raised `FieldError` on every row.
        # `carbon_system` calls `calculate_constants` with these same arguments itself, so
        # the precompute was never anything but a way to get it wrong.
        conditions = (temp_c = gd.temperature[i], sal = gd.salinity[i],
                      pres_bar = gd.pressure[i], PT = gd.phosphate[i],
                      SiT = gd.silicate[i], BT = 415.7, unit = "umol")

        try
            pH[i] = carbon_system(; TA = gd.talk[i], DIC = gd.tco2[i], conditions...).pHtot
        catch
        end
        try
            TA[i] = carbon_system(; pHtot = gd.phtsinsitutp[i], DIC = gd.tco2[i],
                                  conditions...).TA
        catch
        end
        try
            DIC[i] = carbon_system(; pHtot = gd.phtsinsitutp[i], TA = gd.talk[i],
                                   conditions...).DIC
        catch
        end
    end

    return (pH = pH, TA = TA, DIC = DIC)
end


@testset "GLODAPv2" begin
    with_dataset("GLODAPv2", GLODAP_SUBSET, fetch_GLODAP) do path
        gd = load_GLODAP(path)
        @info "GLODAPv2: $(nrow(gd)) rows, $(round(minimum(gd.temperature), digits=1)) to " *
              "$(round(maximum(gd.temperature), digits=1)) °C, salinity " *
              "$(round(minimum(gd.salinity), digits=1)) to $(round(maximum(gd.salinity), digits=1))"

        computed = solve_GLODAP(gd)

        # Tolerances carry headroom over measured agreement (median +0.0029 pH / -1.09 /
        # +1.05 µmol/kg, IQR 0.020 / 7.8 / 7.5), so they catch a regression without being
        # tripped by the noise inherent in field data.
        check_agreement("pH from TA and DIC", gd.phtsinsitutp, computed.pH;
                        max_median = 0.01, max_iqr = 0.05)
        check_agreement("TA from pH and DIC", gd.talk, computed.TA;
                        max_median = 5.0, max_iqr = 20.0)
        check_agreement("DIC from pH and TA", gd.tco2, computed.DIC;
                        max_median = 5.0, max_iqr = 20.0)

        if FIELD_PLOTS
            figures = joinpath(FIGURE_DIR, "GLODAP")
            mkpath(figures)
            savefig(comparison_figure(gd.phtsinsitutp, computed.pH, "pH", "Depth", gd.depth,
                                      lims = (7.4, 8.3), diff_lims = (-0.15, 0.15)),
                    joinpath(figures, "pH_comparison.png"))
            savefig(comparison_figure(gd.talk, computed.TA, "Alk", "Depth", gd.depth,
                                      lims = (1800, 2500), diff_lims = (-60, 60)),
                    joinpath(figures, "TA_comparison.png"))
            savefig(comparison_figure(gd.tco2, computed.DIC, "DIC", "Depth", gd.depth,
                                      lims = (1600, 2400), diff_lims = (-60, 60)),
                    joinpath(figures, "DIC_comparison.png"))
            @info "GLODAPv2: figures written to $figures"
        end
    end
end
