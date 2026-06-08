using DataFrames, CSV, Statistics, Printf, Plots, Measures
using Base.Filesystem
using ProgressMeter
using Downloads
include("../src/CarbonateCalculator.jl")
using .CarbonateCalculator: carbon_system, K_calculator
include("../src/helpers.jl")
using .Helpers
default(
    dpi = 300,
    titlefont  = font(12, "Arial", :darkgray),
    guidefont  = font(10, "Arial", :darkgray),
    tickfont   = font(8,  "Arial", :darkgray),
    legendfont = font(9,  "Arial", :darkgray)
)

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
                 xlims=lims, ylims=lims, xlabel="CODAP Measured", ylabel="Julia predicted")
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]], color=:gray, linestyle=:dash, lw=2)
    
    # Add variable text in the top left
    annotate!(p1, lims[1] + 0.05*(lims[2]-lims[1]), lims[2] - 0.05*(lims[2]-lims[1]),
              text("CODAP_NA_v2021 "*var_name, :left, :top, 12, :darkgray, :bold, "Arial"))

# --- AXIS 2: Measured vs Difference (Residuals) ---
    p2 = scatter(obs_v, diff_v, zcolor=c_v, markerstrokewidth=0, markersize=3,
                 seriesalpha=alpha, legend=false, 
                 colorbar=false,
                 xlims=lims, ylims=diff_lims, xlabel="CODAP Measured", ylabel="measured - predicted",
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


function CODAP_comparison(figdir=".")
    println("\n********************************************")
    println("Generating CODAP-NA v2021 Comparison Plots")
    println("********************************************")

    # 1. Define Local Path and URL
    local_path = joinpath(figdir, "CODAP_NA_v2021.csv")
    url = "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0219960/Data_CODAP/CSV/CODAP_NA_v2021.csv"

    if !isfile(local_path)
        println("CODAP-NA data not found locally. Downloading from NOAA NCEI...")
        try
            # Downloads the file and saves it to your specified local_path
            Downloads.download(url, local_path)
            println("Download complete!")
        catch e
            println("Error downloading the file. Please check your internet connection or the URL.")
            println("Error details: ", e)
            return
        end
    else
        println("Found CODAP-NA file locally. Skipping download.")
    end
    # 2. Load CODAP Data
    println("Importing CODAP-NA Data...")
    # CODAP usually uses standard comma-separation and -999 for missing values
    df = CSV.read(local_path, DataFrame, missingstring=["-999", "-9999", "NaN"])


    # 3. Rename Columns to Standard Internal Names
    rename_dict = Dict(
        "TEMP_Carbonate" => "T",
        "recommended_Salinity_PSS78" => "S",
        "CTDPRES" => "P",
        "DIC" => "DIC",
        "TALK" => "TA",
        "pH_TS_measured" => "pH_obs",
        "Silicate" => "SiT",
        "Phosphate" => "PT"
    )
    
    # Safely rename only the columns that exist
    for (old, new) in rename_dict
        if hasproperty(df, Symbol(old))
            rename!(df, Symbol(old) => Symbol(new))
        end
    end

    function safe_parse(x)
        ismissing(x) && return missing
        parsed = tryparse(Float64, string(x))
        return isnothing(parsed) ? missing : parsed
    end
    
    for col in [:T, :S, :P, :DIC, :TA, :pH_obs, :SiT, :PT]
        if hasproperty(df, col)
            df[!, col] = safe_parse.(df[!, col])
        end
    end

    # Handle missing nutrients: Coalesce turns `missing` into `0.0`
    if hasproperty(df, :SiT)
        df.SiT = coalesce.(df.SiT, 0.0)
    else
        df.SiT = zeros(nrow(df)) # Create empty column if it didn't exist at all
    end

    if hasproperty(df, :PT)
        df.PT = coalesce.(df.PT, 0.0)
    else
        df.PT = zeros(nrow(df))
    end
    
    
    # 4. Clean Data (Require DIC, TA, pH, T, S, P)
    dropmissing!(df, [:DIC, :TA, :pH_obs, :T, :S, :P])
    filter!(row -> row.S > 20 && row.DIC > 1000 && row.TA > 1000, df)
    
    # Convert Pressure from dbar to bar
    df.P .= df.P ./ 10.0

    # 5. Create Figures directory
    fig_folder = joinpath(figdir, "Figures_CODAP")
    isdir(fig_folder) || mkdir(fig_folder)


    # 7. Execute Calculations
    println("Calculating pH, DIC, and TA from different input pairs...")
    
    n_obs    = nrow(df)
    calc_pH  = zeros(Float64, n_obs)
    calc_DIC = zeros(Float64, n_obs)
    calc_TA  = zeros(Float64, n_obs)

    @showprogress for i in 1:nrow(df)
        try
            # --- 1. Calculate constants strictly for THIS row ---
            row_K_results = K_calculator(
                T_in = df.T[i], 
                S_in = df.S[i], 
                P_in = df.P[i],
                K_method = "default"
            )

            # --- 2. Calculate pH (Inputs: TA & DIC) ---
            res_pH = carbon_system(
                TA = df.TA[i], DIC = df.DIC[i],
                T_in = df.T[i], S_in = df.S[i], P_in = df.P[i],
                PT = df.PT[i], SiT = df.SiT[i],
                unit = "umol", Ks = row_K_results.Ks # Pass the extracted tuple!
            )
            calc_pH[i] = res_pH.pHtot

            # --- 3. Calculate DIC (Inputs: TA & pH) ---
            res_DIC = carbon_system(
                TA = df.TA[i], pHtot = df.pH_obs[i],
                T_in = df.T[i], S_in = df.S[i], P_in = df.P[i],
                PT = df.PT[i], SiT = df.SiT[i],
                unit = "umol", Ks = row_K_results.Ks
            )
            calc_DIC[i] = res_DIC.DIC

            # --- 4. Calculate TA (Inputs: DIC & pH) ---
            res_TA = carbon_system(
                DIC = df.DIC[i], pHtot = df.pH_obs[i],
                T_in = df.T[i], S_in = df.S[i], P_in = df.P[i],
                PT = df.PT[i], SiT = df.SiT[i],
                unit = "umol", Ks = row_K_results.Ks
            )
            calc_TA[i] = res_TA.TA

        catch
            calc_pH[i]  = NaN
            calc_DIC[i] = NaN
            calc_TA[i]  = NaN
        end
    end

    # 8. Validation Plots
    println("Making validation plots...")
    
    # pH Plot - Colored by Silicate
    fig_pH = cplot(df.pH_obs, calc_pH, "pH", "Silicate (umol/kg)", df.SiT,
                 lims=(7.2, 8.4), diff_lims=(-0.2, 0.2))
    savefig(fig_pH, joinpath(fig_folder, "CODAP_pH_vs_SiT.png"))

    # DIC Plot - Colored by Silicate
    fig_DIC = cplot(df.DIC, calc_DIC, "DIC", "Silicate (umol/kg)", df.SiT,
                 lims=(1800, 2400), diff_lims=(-50, 50))
    savefig(fig_DIC, joinpath(fig_folder, "CODAP_DIC_vs_SiT.png"))

    # TA Plot - Colored by Silicate
    fig_TA = cplot(df.TA, calc_TA, "TA", "Silicate (umol/kg)", df.SiT,
                 lims=(2000, 2600), diff_lims=(-50, 50))
    savefig(fig_TA, joinpath(fig_folder, "CODAP_TA_vs_SiT.png"))
end