using DataFrames, CSV, Statistics, Printf, Plots, Measures
using Base.Filesystem
using ProgressMeter
using Downloads
include("../src/CarbonateCalculator.jl")
using .CarbonateCalculator: carbon_system, K_calculator
include("../src/helpers.jl")
using .Helpers
using NCDatasets
using DataFrames
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
                 xlims=lims, ylims=lims, xlabel="BGC Argo Measured", ylabel="Julia predicted")
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]], color=:gray, linestyle=:dash, lw=2)
    
    # Add variable text in the top left
    annotate!(p1, lims[1] + 0.05*(lims[2]-lims[1]), lims[2] - 0.05*(lims[2]-lims[1]),
              text("BGC Argo"*var_name, :left, :top, 12, :darkgray, :bold, "Arial"))

# --- AXIS 2: Measured vs Difference (Residuals) ---
    p2 = scatter(obs_v, diff_v, zcolor=c_v, markerstrokewidth=0, markersize=3,
                 seriesalpha=alpha, legend=false, 
                 colorbar=false, 
                 xlims=lims, ylims=diff_lims, xlabel="BGC Argo Measured", ylabel="measured - predicted",
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
    # 1. Generate 200 bins
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


using DataFrames, CSV, HTTP

function BGCArgo_comparison(figdir=".")
    println("\n********************************************")
    println("Generating BGC-Argo Comparison Plots")
    println("********************************************")

    # 1. The ERDDAP API URL
    url = "https://polarwatch.noaa.gov/erddap/tabledap/SOCCOM_BGC_Argo.csv?depth%2Cpressure%2Ctemperature%2Csalinity%2CpH_insitu%2CTALK_LIAR%2CDIC_LIAR&time%3E=2021-10-15T00%3A00%3A00Z&time%3C=2021-10-22T14%3A47%3A00Z"

    local df 

    println("Requesting data from ERDDAP...")
    
    try
        resp = HTTP.get(url)
        df_raw = CSV.read(IOBuffer(resp.body), DataFrame, header=1, skipto=3)
        
        if nrow(df_raw) == 0
            println("Server returned 0 rows. Try widening the time range!")
            return
        end

        # Map the columns and add blank nutrient columns for the calculator
        df = DataFrame(
            depth  = Float64.(df_raw[:, 1]),
            P      = Float64.(df_raw[:, 2]),
            T      = Float64.(df_raw[:, 3]),
            S      = Float64.(df_raw[:, 4]),
            pH_obs = Float64.(df_raw[:, 5]),
            TA     = Float64.(df_raw[:, 6]),
            DIC    = Float64.(df_raw[:, 7]),
            SiT    = zeros(Float64, nrow(df_raw)), # Placeholder for Silicate
            PT     = zeros(Float64, nrow(df_raw))  # Placeholder for Phosphate
        )
        
    catch e
        println("An unexpected error occurred during download: ", e)
        return
    end

    println("Success! Captured $(nrow(df)) data points.")
    
    # 2. Clean the Data
    # Standardize missing values (Argo often uses -999 for bad data)
    df.T .= [ (x < -5 || x > 40) ? missing : x for x in df.T ]
    df.S .= [ (x < 0 || x > 50) ? missing : x for x in df.S ]
    df.P .= [ (x < 0 || x > 6000) ? missing : x for x in df.P ]

    dropmissing!(df, [:P, :T, :S, :pH_obs, :DIC, :TA])
    
    # Convert Pressure from dbar to bar 
    df.P .= df.P ./ 10.0

    println("Cleaned data down to $(nrow(df)) valid points. Ready for calculations!")

    # 3. Create Figures directory
    fig_folder = joinpath(figdir, "Figures_BGCArgo")
    isdir(fig_folder) || mkdir(fig_folder)

    # 5. Execute Calculations
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

            # --- 2. Calculate pH ---
            res_pH = carbon_system(
                TA = df.TA[i], DIC = df.DIC[i],
                T_in = df.T[i], S_in = df.S[i], P_in = df.P[i],
                PT = df.PT[i], SiT = df.SiT[i],
                unit = "umol", Ks = row_K_results.Ks # Pass the extracted tuple!
            )
            calc_pH[i] = res_pH.pHtot

            # --- 3. Calculate DIC ---
            res_DIC = carbon_system(
                TA = df.TA[i], pHtot = df.pH_obs[i],
                T_in = df.T[i], S_in = df.S[i], P_in = df.P[i],
                PT = df.PT[i], SiT = df.SiT[i],
                unit = "umol", Ks = row_K_results.Ks
            )
            calc_DIC[i] = res_DIC.DIC

            # --- 4. Calculate TA ---
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

    # 6. Validation Plots
    println("Making validation plots...")
    
    fig_pH = cplot(df.pH_obs, calc_pH, "pH", "Depth (m)", df.depth,
                 lims=(7.4, 8.2), diff_lims=(-0.1, 0.1))
    savefig(fig_pH, joinpath(fig_folder, "BGCArgo_pH_Comparison.png"))

    fig_DIC = cplot(df.DIC, calc_DIC, "DIC", "Depth (m)", df.depth,
                 lims=(1900, 2400), diff_lims=(-30, 30))
    savefig(fig_DIC, joinpath(fig_folder, "BGCArgo_DIC_Comparison.png"))

    fig_TA = cplot(df.TA, calc_TA, "TA", "Depth (m)", df.depth,
                 lims=(2200, 2500), diff_lims=(-30, 30))
    savefig(fig_TA, joinpath(fig_folder, "BGCArgo_TA_Comparison.png"))
    
    println("Finished! Check the Figures_BGCArgo folder for your deep-water plots.")
end