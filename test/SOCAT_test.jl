using DataFrames, CSV, Statistics, Printf, Plots, Measures, Downloads
using Base.Filesystem
using ProgressMeter
using Statistics
include("../src/Calculator.jl")
using .Calculator: carbon_system, K_calculator
include("../src/helpers.jl")
using .Helpers
default(
    dpi = 300,
    titlefont  = font(12, "Arial", :darkgray),
    guidefont  = font(10, "Arial", :darkgray),
    tickfont   = font(8,  "Arial", :darkgray),
    legendfont = font(9,  "Arial", :darkgray)
)

# ------------------------------------------------------------------- #
# Helper: Simple Alkalinity Proxy for SOCAT
# ------------------------------------------------------------------- #
# Since SOCAT doesn't have TA, we estimate it from Salinity (S)
# Using a global average (Millero et al. 2006 style)
function estimate_TA(S, T)
    return 2300.0 * (S / 35.0)
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
                 xlims=lims, ylims=lims, xlabel="SOCAT Measured", ylabel="Julia predicted")
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]], color=:gray, linestyle=:dash, lw=2)
    
    # Add variable text in the top left
    annotate!(p1, lims[1] + 0.05*(lims[2]-lims[1]), lims[2] - 0.05*(lims[2]-lims[1]),
              text("SOCAT "*var_name, :left, :top, 12, :darkgray, :bold, "Arial"))

# --- AXIS 2: Measured vs Difference (Residuals) ---
    p2 = scatter(obs_v, diff_v, zcolor=c_v, markerstrokewidth=0, markersize=3,
                 seriesalpha=alpha, legend=false, 
                 colorbar=false,
                 xlims=lims, ylims=diff_lims, xlabel="SOCAT Measured", ylabel="measured - predicted",
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


function SOCAT_comparison(figdir=".")
    socat_url = "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0304549/SOCATv2025_NorthPacific.tsv"
    local_path = joinpath(figdir, "SOCATv2025_NorthPacific.tsv")

    if !isfile(local_path)
        println("Fetching North Pacific SOCAT from NCEI...")
        prog = Progress(100, 1, "Downloading: ")

        Downloads.download(socat_url, local_path; timeout=600.0, progress = (total, now) -> begin
            if total > 0
                update!(prog, round(Int, (now / total) * 100))
            end
        end)
    else
        println("Found file locally.")
    end


    println("Streaming SOCAT Data (Systematic Sampling)...")

    data_start = 0
    open(local_path) do f
        for (i, line) in enumerate(eachline(f))
            if occursin("Expocode", line) && (occursin("fCO2", line) || occursin("temp", line))
                data_start = i + 1 # The actual numbers start on the line AFTER this
                break
            end
        end
    end

    if data_start == 0
        error("Could not find the technical data header!")
    end

    println("Skipping metadata. Data starts on line $data_start. Starting stream...")

    df = DataFrame()
    row_count = 0
    sampled_count = 0

    # Skip the top of the file (metadata)
    rows = CSV.Rows(local_path,
                    delim='\t',
                    comment="#",
                    header=false,
                    skipto=data_start)

    for row in rows
        row_count += 1
        if row_count % 1000 == 0
            push!(df, row, cols=:union)
            sampled_count += 1
        end
    end


    # Col 14 = Salinity, Col 15 = Temp, Col 23 = fCO2
    rename!(df, :Column14 => :S, :Column15 => :T, :Column23 => :fCO2)

    # Clean up the data types to ensure they are numbers, not strings
    df.S = parse.(Float64, string.(df.S))
    df.T = parse.(Float64, string.(df.T))
    df.fCO2 = parse.(Float64, string.(df.fCO2))

    # Drop missing values and perform basic quality control
    dropmissing!(df, [:fCO2, :T, :S])
    filter!(row -> row.fCO2 > 0 && row.S > 5, df)

    # Safety check before calculating!
    if nrow(df) == 0
        error("The dataframe is empty after filtering! The columns might still be misaligned.")
    else
        println("Successfully loaded and filtered $(nrow(df)) valid observations.")
    end


    n_obs = nrow(df)
    calc_fCO2_out = zeros(Float64, n_obs)
    calc_pH       = zeros(Float64, n_obs)
    calc_DIC      = zeros(Float64, n_obs)


    println("Calculating System from fCO2 and Estimated TA...")

    @showprogress for i in 1:nrow(df)
        try
            # --- 1. Calculate constants strictly for THIS row ---
            # (Assuming P=0 since SOCAT is strictly surface underway data)
            row_K_results = K_calculator(
                T_in = df.T[i], 
                S_in = df.S[i], 
                P_in = 0.0,
                K_method = "default"
            )

            # --- 2. Calculate System ---
            # Use fCO2 and estimated TA as our two knowns
            res = carbon_system(
                fCO₂ = df.fCO2[i],
                TA = estimate_TA(df.S[i], df.T[i]),
                T_in = df.T[i],
                S_in = df.S[i],
                P_in = 0.0,
                unit = "umol",
                Ks = row_K_results.Ks
            )
            
            calc_fCO2_out[i] = res.fCO₂
            calc_pH[i] = res.pHtot
            calc_DIC[i] = res.DIC
            
        catch e
            calc_fCO2_out[i] = NaN
            calc_pH[i] = NaN
            calc_DIC[i] = NaN
        end
    end

    println("Making validation plots...")
    fig1 = cplot(df.fCO2, calc_fCO2_out, "fCO2", "Temp", df.T, 
                 lims=(200, 600), diff_lims=(-0.5, 0.5))
    
    savefig(fig1, joinpath(figdir, "North_Pacific_fCO2_RoundTrip.png"))
    
    println("Finished! Check the SOCAT Figures folder.")
end

