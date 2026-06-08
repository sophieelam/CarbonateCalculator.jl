using Downloads
using ZipFile
using DataFrames
using CSV
using ProgressMeter
using Base.Filesystem
using Plots
using Statistics
using Printf
include("../src/CarbonateCalculator.jl")
using .CarbonateCalculator: carbon_system, boron_system, boron_isotopes, whole_system
using PythonCall
const np = pyimport("numpy")
const kgen = pyimport("kgen")
include("../src/helpers.jl")
using .Helpers
default(
    dpi = 300,
    titlefont  = font(12, "Arial", :darkgray),
    guidefont  = font(10, "Arial", :darkgray), # Axis labels
    tickfont   = font(8,  "Arial", :darkgray), # Numbers on the axes
    legendfont = font(9,  "Arial", :darkgray)
)


"""
Downloads, extracts, and cleans GLODAPv2 data.
"""
function get_GLODAP(; path=".", leave_zip=true)
    zip_path = joinpath(path, "GLODAPv2 Merged Master File.csv.zip")
    
    if !isfile(zip_path)
        println("Fetching GLODAPv2 Data...")

        GLODAP_urls = [
            "https://www.glodap.info/glodap_files/v2.2020/GLODAPv2.2020_Merged_Master_File.csv.zip"
            # "https://www.glodap.info/glodap_files/v2.2023/GLODAPv2.2023_Merged_Master_File.csv.zip"
        ]
        
        downloaded = false
        for url in GLODAP_urls
            if downloaded
                break
            end
            
            try
                prog = Progress(100, 1, "Downloading GLODAPv2: ")
                
                Downloads.download(url, zip_path; progress = (total, now) -> begin
                    if total > 0
                        update!(prog, round(Int, (now / total) * 100))
                    end
                end)
                
                finish!(prog)
                downloaded = true
            catch e
                println("Failed to download from $url. Trying next...")
            end
        end
    else
        println("Found GLODAPv2 Data...")
    end

    println("Reading data...")
    z = ZipFile.Reader(zip_path)

    target_filename = "GLODAPv2.2020_Merged_Master_File.csv"
    file_in_zip = first(filter(f -> f.name == target_filename, z.files))
    
    # Read data directly into DataFrame; convert -9999 to `missing`!
    gd = CSV.read(read(file_in_zip), DataFrame, missingstrings=["-9999", "-9999.0", "-999.9"])
    close(z)

    println("Selecting 'good' (flag == 2) data...")
    # Isolate good data only (flag = 2). 
    gd.phtsinsitutp = ifelse.(gd.phtsinsitutpf .== 2, gd.phtsinsitutp, missing)
    gd.tco2         = ifelse.(gd.tco2f .== 2, gd.tco2, missing)
    gd.talk         = ifelse.(gd.talkf .== 2, gd.talk, missing)
    gd.salinity     = ifelse.(gd.salinityf .== 2, gd.salinity, missing)
    gd.phosphate    = ifelse.(gd.phosphatef .== 2, gd.phosphate, missing)
    gd.silicate     = ifelse.(gd.silicatef .== 2, gd.silicate, missing)

    # Drop rows where any of the critical columns are missing
    cols_to_check = [:phtsinsitutp, :tco2, :talk, :temperature, :salinity, :pressure, :silicate, :phosphate]
    dropmissing!(gd, cols_to_check)

    println("Saving data subset...")
    
    # Select columns to keep
    cols_to_keep = [
        :phts25p0, :phtsinsitutp, :tco2, :talk, :temperature, :salinity,
        :cruise, :station, :cast, :year, :month, :day, :hour,
        :latitude, :longitude, :bottomdepth, :maxsampdepth, :bottle,
        :pressure, :depth, :theta, :silicate, :phosphate
    ]
    gds = select(gd, cols_to_keep)

    # Save to CSV
    subset_path = joinpath(path, "GLODAPv2_pH_DIC_ALK_subset.csv")
    CSV.write(subset_path, gds)

    if !leave_zip
        rm(zip_path)
    end

    return nothing
end


function cplot(obs, pred, var_name::String, cvar_name::String, c_data; alpha=0.4, pclims=[0.05, 0.9995], lims=nothing, diff_lims=nothing, hist_xlims=nothing)
    
    # Calculate difference and remove NaNs
    diff = obs .- pred
    valid = .!isnan.(obs) .& .!isnan.(pred)
    obs_v, pred_v, diff_v, c_v = obs[valid], pred[valid], diff[valid], c_data[valid]

    # Limits for Measured vs Predicted (Plot 1)
    if lims === nothing
        ad = vcat(obs_v, pred_v)
        mn, mx = quantile(ad, pclims)
        pad = 0.1 * (mx - mn)
        lims = (mn - pad, mx + pad)
    end

    # Limits for Residuals (Plots 2 & 3)
    if diff_lims === nothing
        diff_mn, diff_mx = quantile(diff_v, pclims)
        diff_pad = 0.15 * (diff_mx - diff_mn)
        diff_lims = (diff_mn - diff_pad, diff_mx + diff_pad)
    end

    # --- AXIS 1: Measured vs Predicted ---
    p1 = scatter(obs_v, pred_v, zcolor=c_v, markerstrokewidth=0, markersize=3,
                 seriesalpha=alpha, legend=false, colorbar=false, 
                 xlims=lims, ylims=lims, xlabel="GLODAPv2 Measured", ylabel="Julia predicted")
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]], color=:gray, linestyle=:dash, lw=2)
    
    # Add variable text in the top left
    annotate!(p1, lims[1] + 0.05*(lims[2]-lims[1]), lims[2] - 0.05*(lims[2]-lims[1]),
              text("GLODAPv2 "*var_name, :left, :top, 12, :darkgray, :bold, "Arial"))

# --- AXIS 2: Measured vs Difference (Residuals) ---
    p2 = scatter(obs_v, diff_v, zcolor=c_v, markerstrokewidth=0, markersize=3,
                 seriesalpha=alpha, legend=false, 
                 colorbar=false,
                 xlims=lims, ylims=diff_lims, xlabel="GLODAPv2 Measured", ylabel="measured - predicted",
                 right_margin=2Plots.mm)
    hline!(p2, [0], color=:gray, linestyle=:dash, lw=2)

    # Calculate Stats for Annotations
    med_diff = median(diff_v)
    lim95 = quantile(diff_v, [0.025, 0.975])
    
    fmt(val) = abs(val) < 0.01 ? @sprintf("%.1e", val) : @sprintf("%.2f", val)
    stat_text = "Median Offset: $(fmt(med_diff))\n95% Limits: $(fmt(lim95[1] - med_diff)) / +$(fmt(lim95[2] - med_diff))"

    annotate!(p2, lims[1] + 0.03*(lims[2]-lims[1]), diff_lims[2] - 0.03*(diff_lims[2]-diff_lims[1]),
              text(stat_text, :left, :top, 9, :black, "Arial"))

    # --- AXIS 3: Histogram ---
    bin_edges = range(diff_lims[1], diff_lims[2], length=200)

    p3 = histogram(diff_v, orientation=:horizontal, bins=bin_edges, color=:gray,
                   legend=false, ylims=diff_lims, xlabel="n", yticks=false,
                   left_margin=0Plots.mm)

    hline!(p3, [med_diff], color=:red, linestyle=:dash, lw=2)
    hspan!(p3, [lim95[1], lim95[2]], color=:red, alpha=0.2)
    hline!(p3, [0], color=:gray, linestyle=:dash, lw=2)

    scatter!(p3, [0], [0], zcolor=[c_v[1]], clims=(minimum(c_v), maximum(c_v)),
             markeralpha=0, markersize=0, label="",
             colorbar=true, colorbar_title=cvar_name)

    # Enforce y-limits to prevent Plots.jl cropping
    ylims!(p2, diff_lims)
    ylims!(p3, diff_lims)

    # --- Combine into Layout ---
    l = @layout [a{0.35w} b{0.40w} c{0.25w}]
    fig = plot(p1, p2, p3, layout=l, size=(1200, 450), margin=6Plots.mm)
    return fig
end



function carbon_systemGLODAPv2_comparison(figdir=".")
    println("\n********************************************")
    println("Generating GLODAPv2 Comparison Plots")
    println("********************************************")

    filepath = joinpath(figdir, "GLODAPv2_pH_DIC_ALK_subset.csv")
    
    if !isfile(filepath)
        println("Warning: Data file not found at $filepath")
        println("Triggering download and extraction process now...")
        get_GLODAP(path=figdir)
    end


    println("Importing GLODAPv2 Data...")
    gd = CSV.read(filepath, DataFrame)

    # Drop rows missing any of the critical variables
    cols_to_check = [:phtsinsitutp, :temperature, :salinity, :tco2, :talk, :pressure, :phosphate, :silicate]
    dropmissing!(gd, cols_to_check)

    filter!(row -> row.salinity > 20.0 && row.tco2 > 1000.0 && row.talk > 1000.0, gd)

    # Convert pressure to bar
    gd.pressure .= gd.pressure ./ 10.0

    # Exclude weird cruise 270 data
    filter!(row -> row.cruise != 270, gd)

    fig_folder = joinpath(figdir, "Figures")
    isdir(fig_folder) || mkdir(fig_folder)


    println("Pre-calculating Constants...")


    # ================================================================= #
    # --- CALCULATION 1: pH from DIC and TA ---
    # ================================================================= #
    println("Calculating pH from DIC and TA...")
    cpH_pHtot = zeros(nrow(gd))
    for i in 1:nrow(gd)
        try
            # 1. Calculate the constants for just THIS row
            row_K_results = CarbonateCalculator.K_calculator(
                T_in = gd.temperature[i],
                S_in = gd.salinity[i],
                P_in = gd.pressure[i],
                K_method = "default"
            )

            # 2. Pass the purely scalar inputs (and the inner .Ks bundle) to the system
            res = carbon_system(
                TA = gd.talk[i],
                DIC = gd.tco2[i],
                T_in = gd.temperature[i],
                S_in = gd.salinity[i],
                P_in = gd.pressure[i],
                PT = gd.phosphate[i],
                SiT = gd.silicate[i],
                BT = 415.7, 
                unit = "umol", 
                Ks = row_K_results.Ks
            )
            cpH_pHtot[i] = res.pHtot
            
        catch e
            cpH_pHtot[i] = NaN
        end
    end
    
    println("Failed rows: ", count(isnan, cpH_pHtot))
    println("  Making plots...")


    fig1 = cplot(gd.phtsinsitutp, cpH_pHtot, "pH", "Depth", gd.depth, lims=(7.4, 8.3), diff_lims=(-0.15, 0.15), hist_xlims=(0, 10000))
    savefig(fig1, joinpath(fig_folder, "pH_comparison.png"))

    # --- CALCULATION 2: TA from pH and DIC ---
    println("Calculating TA from pH and DIC...")
    cTA_TA = zeros(nrow(gd))
    
    for i in 1:nrow(gd)
        try
            # 1. Calculate the constants for just THIS row
            row_K_results = CarbonateCalculator.K_calculator(
                T_in = gd.temperature[i],
                S_in = gd.salinity[i],
                P_in = gd.pressure[i],
                K_method = "default"
            )

            # 2. Pass the purely scalar inputs (and the inner .Ks bundle) to the system
            res = carbon_system(
                pHtot = gd.phtsinsitutp[i],
                DIC = gd.tco2[i],
                T_in = gd.temperature[i],
                S_in = gd.salinity[i],
                P_in = gd.pressure[i],
                PT = gd.phosphate[i],
                SiT = gd.silicate[i],
                BT = 415.7, 
                unit = "umol", 
                Ks = row_K_results.Ks
            )
            cTA_TA[i] = res.TA
            
        catch e
            cTA_TA[i] = NaN
        end
    end
    
    println("  Making plots...")
    fig2 = cplot(gd.talk, cTA_TA, "Alk", "Depth", gd.depth, lims=(1800, 2500), diff_lims=(-60, 60), hist_xlims=(0, 10000))
    savefig(fig2, joinpath(fig_folder, "TA_comparison.png"))


    # --- CALCULATION 3: DIC from pH and TA ---
    println("Calculating DIC from pH and TA...")
    cDIC_DIC = zeros(nrow(gd))
    
    for i in 1:nrow(gd)
        try
            # 1. Calculate the constants for just THIS row
            row_K_results = CarbonateCalculator.K_calculator(
                T_in = gd.temperature[i],
                S_in = gd.salinity[i],
                P_in = gd.pressure[i],
                K_method = "default"
            )

            # 2. Pass the purely scalar inputs (and the inner .Ks bundle) to the system
            res = carbon_system(
                pHtot = gd.phtsinsitutp[i],
                TA = gd.talk[i],
                T_in = gd.temperature[i],
                S_in = gd.salinity[i],
                P_in = gd.pressure[i],
                PT = gd.phosphate[i],
                SiT = gd.silicate[i],
                BT = 415.7, 
                unit = "umol", 
                Ks = row_K_results.Ks # Extracted right from the row calculation!
            )
            cDIC_DIC[i] = res.DIC
            
        catch e
            cDIC_DIC[i] = NaN
        end
    end
    
    println("  Making plots...")
    fig3 = cplot(gd.tco2, cDIC_DIC, "DIC", "Depth", gd.depth, lims=(1600, 2400), diff_lims=(-60, 60), hist_xlims=(0, 10000))
    savefig(fig3, joinpath(fig_folder, "DIC_comparison.png"))

end


if abspath(PROGRAM_FILE) == @__FILE__
    carbon_systemGLODAPv2_comparison()
end