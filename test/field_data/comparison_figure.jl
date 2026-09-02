# Shared three-panel comparison figure for block 3, loaded only when CC_PLOTS=true and
# Plots is installed. All four datasets plot the same thing — measured against computed,
# the residuals, and their distribution — differing only in what colours the points.

"""
Three-panel comparison: measured against computed, the residuals against measured, and
the residual distribution, coloured by `colour_data`.
"""
function comparison_figure(measured, computed, var_name, colour_name, colour_data;
                           dataset = "measured", alpha = 0.4, lims = nothing,
                           diff_lims = nothing)
    residual = measured .- computed
    valid = .!isnan.(measured) .& .!isnan.(computed)
    measured_v, computed_v = measured[valid], computed[valid]
    residual_v, colour_v = residual[valid], colour_data[valid]

    if lims === nothing
        low, high = quantile(vcat(measured_v, computed_v), [0.05, 0.9995])
        pad = 0.1 * (high - low)
        lims = (low - pad, high + pad)
    end
    if diff_lims === nothing
        low, high = quantile(residual_v, [0.05, 0.9995])
        pad = 0.15 * (high - low)
        diff_lims = (low - pad, high + pad)
    end

    p1 = scatter(measured_v, computed_v, zcolor = colour_v, markerstrokewidth = 0,
                 markersize = 3, seriesalpha = alpha, legend = false, colorbar = false,
                 xlims = lims, ylims = lims,
                 xlabel = "$dataset measured", ylabel = "CarbonateCalculator.jl")
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]],
          color = :gray, linestyle = :dash, lw = 2)
    annotate!(p1, lims[1] + 0.05 * (lims[2] - lims[1]),
              lims[2] - 0.05 * (lims[2] - lims[1]),
              text("$dataset $var_name", :left, :top, 12, :darkgray, :bold))

    p2 = scatter(measured_v, residual_v, zcolor = colour_v, markerstrokewidth = 0,
                 markersize = 3, seriesalpha = alpha, legend = false, colorbar = false,
                 xlims = lims, ylims = diff_lims, xlabel = "$dataset measured",
                 ylabel = "measured - computed", right_margin = 2Plots.mm)
    hline!(p2, [0], color = :gray, linestyle = :dash, lw = 2)

    median_residual = median(residual_v)
    limits95 = quantile(residual_v, [0.025, 0.975])
    fmt(v) = abs(v) < 0.01 ? @sprintf("%.1e", v) : @sprintf("%.2f", v)
    annotate!(p2, lims[1] + 0.03 * (lims[2] - lims[1]),
              diff_lims[2] - 0.03 * (diff_lims[2] - diff_lims[1]),
              text("Median offset: $(fmt(median_residual))\n95% limits: " *
                   "$(fmt(limits95[1] - median_residual)) / " *
                   "+$(fmt(limits95[2] - median_residual))", :left, :top, 9, :black))

    p3 = histogram(residual_v, orientation = :horizontal,
                   bins = range(diff_lims[1], diff_lims[2], length = 200),
                   color = :gray, legend = false, ylims = diff_lims,
                   xlabel = "n", yticks = false, left_margin = 0Plots.mm)
    hline!(p3, [median_residual], color = :red, linestyle = :dash, lw = 2)
    hspan!(p3, [limits95[1], limits95[2]], color = :red, alpha = 0.2)
    hline!(p3, [0], color = :gray, linestyle = :dash, lw = 2)
    scatter!(p3, [0], [0], zcolor = [colour_v[1]],
             clims = (minimum(colour_v), maximum(colour_v)), markeralpha = 0,
             markersize = 0, label = "", colorbar = true,
             colorbar_title = colour_name)

    return plot(p1, p2, p3, layout = Plots.grid(1, 3, widths = [0.35, 0.40, 0.25]),
                size = (1200, 450), margin = 6Plots.mm)
end
