# --- Publication figure settings --------------------------------------------
# Target journal requirements:
#   * vector PDF, 300-600 dpi (a vector PDF has effectively unlimited resolution)
#   * single-column width 8.5 cm when possible; absolute max 18 cm wide x 22 cm tall
#   * all text between 6 and 10 pt at publication size
#
# The GR backend writes PDF with 1 `size` unit = 1 PostScript point (1 pt = 1/72"),
# and font sizes are absolute points, independent of the figure size (`dpi` only
# affects raster output). So figure sizes are given directly in points, derived
# from the physical dimensions below, and every font size is kept within 6-10 pt.
const _PT_PER_CM = 72 / 2.54       # ≈ 28.35 pt per cm
_cm(x) = round(Int, x * _PT_PER_CM)
const _MAX_W = _cm(18)             # 510 pt — absolute maximum width  (18 cm)
const _MAX_H = _cm(22)             # 624 pt — absolute maximum height (22 cm)
const _COL_W = _cm(8.5)            # 241 pt — single-column width     (8.5 cm)

Plots.default(
    guidefontsize      = 8,
    tickfontsize       = 7,
    legendfontsize     = 6,
    titlefontsize      = 9,
    annotationfontsize = 8,
    dpi                = 600,
)

function _plot_early_warning_signals(timestamp, p, s, variance_ts, autocorr_ts, time_horizons, times,
                                     threshold_idx, kendall_tau)

    pd_label = map(idx -> "Time horizon = $(time_horizons[idx]). Kendall-τ = $(kendall_tau[Symbol(:s_, idx)])",
                   1:length(time_horizons)
                  )
    pd_label = reshape(pd_label, (1, size(s, :time_horizon)))
    xticks = 0:25:length(times)
    xlims = (0, times[end]+1)

    plt1 = plot(collect(times), p, label=false, ylabel="Amount of Phosphorus")
    vline!([times[threshold_idx]], label="Threshold", color="black", lw=2, xticks=xticks, xlims=xlims)
    label1 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(a)"))

    plt2 = plot(collect(times), s, label=pd_label, ylabel="Pathway Diversity", legend=:bottomleft)
    vline!([times[threshold_idx]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims)
    label2 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(b)"))

    new_time = times[1 + length(times) - length(variance_ts):end]
    plt3 = plot(collect(new_time), variance_ts, label=label="Kendall-τ=$(kendall_tau[:var])", ylabel="Variance")
    vline!([times[threshold_idx]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims)
    label3 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(c)"))

    new_time = times[1 + length(times) - length(autocorr_ts):end]
    plt4 = plot(collect(new_time), autocorr_ts, label="Kendall-τ=$(kendall_tau[:autocorr])", ylabel="Autocorrelation")
    vline!([times[threshold_idx]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims, xlabel = "Time (year)")
    label4 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(d)"))

    plot(label1, label2, label3, label4, plt1, plt2, plt3, plt4, layout=@layout([grid(4,1){0.005w} grid(4,1)]),
         size=(_MAX_W, _cm(21)), top_margin = 7Plots.mm)
    savefig("../output/$(timestamp)_early_warning_signal.pdf")
end

function _plot_distance_threshold(s, s_diff, P0_options, time_horizons, timestamp, peaks_idx, one_plot)
    if size(s, :time_horizon) == 1
        label = false
    else
        label = reshape(map(t_horizon -> "Time horizon = $(t_horizon)", time_horizons), (1,size(s, :time_horizon)))
    end

    peak_values = [s_diff[P0=peak_idx, time_horizon=idx] for (idx,peak_idx) in enumerate(peaks_idx)]
    xlims = (0, P0_options[end] + 0.1)

    if one_plot
        plt1 = plot(P0_options, s, label=label, legend=:topright, xlabel="Initial Amount of Phosphorus")
        vspan!([0.64, 1.86], linecolor = :grey, fillcolor = :grey, alpha = 0.35, label = false)
        vline!([1.25], label="Median threshold", color="black", ylabel="Pathway Diversity", xlims=xlims)
        annotate!([(0.35,0.16,"Oligotrophic \nbasin of \nattraction"), (2.5,0.16,"Eutrophic \nbasin of attraction"),
                   (0.95,0.14,"Threshold\nregion")])
        label1 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(a)"))

        plt2 = plot(P0_options[2:end], s_diff, label=label, xlabel="Initial Amount of Phosphorus", legend=:topright)
        scatter!(P0_options[2:end][peaks_idx], peak_values, label=false, markerstrokewidth=0, color=[1,2,3,4])
        vspan!([0.64, 1.86], linecolor = :grey, fillcolor = :grey, alpha = 0.35, label = false)
        vline!([1.25], label="Median threshold", color="black", ylabel="Pathway Diversity Derivative", xlims=xlims)
        label2 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(b)"))

        plot(label1, label2, plt1, plt2, layout=@layout([grid(2,1){0.03w} grid(2,1)]), size=(_MAX_W, _cm(20)))
        savefig("../output/$(timestamp)_distance_threshold.pdf")
    else
        plot(P0_options, s, label=label, legend=:topright, xlabel="Amount of Phosphorus")
        vspan!([0.64, 1.86], linecolor = :grey, fillcolor = :grey, alpha = 0.35, label = false)
        vline!([1.25], label="Median threshold", color="black", ylabel="Pathway Diversity",
               left_margin = 10Plots.mm, size=(_cm(14), _cm(9)), xlims=xlims)
        annotate!([(0.3,0.16,"Oligotrophic \nbasin of \nattraction"), (2.5,0.16,"Eutrophic\nbasin of attraction"),
                   (1.0,0.14,"Threshold\nregion")])
        savefig("../output/$(timestamp)_distance_threshold.pdf")

        plot(P0_options[2:end], s_diff, label=label, legend=:topleft, xlabel="Amount of Phosphorus")
        vspan!([0.64, 1.86], linecolor = :grey, fillcolor = :grey, alpha = 0.35, label = false)
        scatter!(P0_options[2:end][peaks_idx], peak_values, label=false, markerstrokewidth=0, color=[1,2,3,4])
        vline!([1.25], label="Median threshold", color="black", ylabel="Pathway Diversity Derivative",
               left_margin = 10Plots.mm, size=(_cm(14), _cm(9)), xlims=xlims)
        savefig("../output/$(timestamp)_distance_threshold_derivative.pdf")
    end
end

function _plot_sensitivity(s, P0_options, timestamp, scenarios, relative = false)
    if relative
        first_idx = 2
        plot(P0_options, s[type=first_idx], label="$(scenarios[first_idx][:name])", legend=:topright,
             ylabel="Relative Pathway Diversity", left_margin = 10Plots.mm, size=(_MAX_W, _cm(12)),
             alpha=0.9, xlabel="Initial Amount of Phosphorus", ylims=[0.67, 1.3])
        vspan!([0.64, 1.86], linecolor = :grey, fillcolor = :grey, alpha = 0.35, label = false)
        annotate!([(0.3,0.85,"Oligotrophic\n basin of\nattraction"), (2.5,0.75,"Eutrophic\n basin of attraction"),
                   (1.0,0.75,"Threshold\nregion")])
        hline!([1.0], label=false, color="black")
    else
        first_idx = 1
        plot(P0_options, s[type=first_idx], label="$(scenarios[first_idx][:name])", legend=:topright,
             ylabel="Pathway Diversity", left_margin = 10Plots.mm, size=(_MAX_W, _cm(12)), color="black", lw=1.5,
             alpha=0.9, xlabel="Initial Amount of Phosphorus")
        vspan!([0.64, 1.86], linecolor = :grey, fillcolor = :grey, alpha = 0.35, label = false)
        annotate!([(0.3,0.12,"Oligotrophic basin\n of attraction"), (2.5,0.16,"Eutrophic\n basin of attraction"),
                   (1.0,0.12,"Threshold\nregion")])
    end

    for (idx_scenario, scenario) in enumerate(scenarios[first_idx+1:end])
        plot!(P0_options, s[type=idx_scenario+first_idx], label="$(scenario[:name])",
              color=idx_scenario+first_idx-1, alpha=0.9)
    end
    vline!([1.25], label="Median threshold", color="black")
    savefig("../output/$(timestamp)_sensitivity.pdf")
end

function _plot_states_distribution(P0_options, n_decision, timestamp)
    plot_array = Plots.Plot[]
    labels = reshape(map(decision -> "Decision $(decision)", 1:n_decision), 1, n_decision)
    letters = ["(a) Very low", "(b) Low", "(c) High", "(d) Very high"]

    for (idx, P0) in enumerate(P0_options)
        states = readdlm("../output/$(timestamp)_state_distribution_$(idx).csv", ',')
        states = [_clean_vector(row) for row in eachrow(states)]
        plt1 = violin(labels, states[2:end], ylim=(0.0, 2.5), ylabel="State", color=idx, label=false,
                      linewidth=0.5, title="$(letters[idx]) initial condition (x=$(round(P0, digits=2)))")
        push!(plot_array,plt1)
    end
    plot(plot_array..., layout=(length(P0_options), 1),
         size=(_cm(17), min(_MAX_H, _cm(5)*length(P0_options))))
    savefig("../output/$(timestamp)_state_distribution.pdf")
end

function _plot_early_warning_residuals(timestamp, p, residuals, times, threshold_idx)
    xticks = 0:5:length(times)
    xlims = (0, times[end]+1)

    plt1 = plot(collect(times), p, label=false, ylabel="Amount of Phosphorus")
    vline!([times[threshold_idx]], label="Threshold", color="black", lw=2, xticks=xticks, xlims=xlims)
    label1 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(a)"))

    plt2 = plot(collect(times), residuals, label=false, ylabel="Residual", xlabel="Time (year)")
    vline!([times[threshold_idx]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims)
    hline!([0.0], label=false, color="black", lw=1)
    label2 = plot(grid = false, showaxis = false, annotation=(0.1,0.5,"(b)"))

    plot(label1, label2, plt1, plt2, layout=@layout([grid(2,1){0.01w} grid(2,1)]),
         size=(_MAX_W, _cm(10)), bottom_margin = 5Plots.mm)
    savefig("../output/$(timestamp)_early_warning_signal_residuals.pdf")
end

function _plot_decision_scales(s, time_horizons, decision_steps, timestamp)
    label = map(decision_step -> "Time horizon = $(decision_step)", time_horizons)
    label = reshape(label, (1,length(time_horizons)))

    plot(decision_steps, s, label = label, legend=:topright, size=(_cm(12), _cm(8)),
         ylabel = "Pathway diversity", xlabel = "Decision step", left_margin = 10Plots.mm)
    savefig("../output/$(timestamp)_decision_scales.pdf")
end


function _plot_scaling(s, P0_options, number_options, timestamp)
    label = map(P0 -> "Initial condition = $(P0)", P0_options)
    label = reshape(label, (1,length(P0_options)))

    plot(number_options, s, label=label, left_margin = 5Plots.mm, legend=:outerbottomright,
         size=(_cm(14), _cm(9)), ylabel = "Pathway diversity", xlabel = "Maximum number of options")
    savefig("../output/$(timestamp)_scaling.pdf")
end

function _plot_number_options_simulation(P0_options, influx_options, n_options)
    heatmap(influx_options, P0_options, n_options, label=false, ylabel="Amount of Phosphorus", xlabel="Influx",
            colorbar_title="Future number of options", left_margin = 5Plots.mm, size=(_cm(12), _cm(9)))
    savefig("../output/number_options_simulation.pdf")
end

function _plot_range_states(P0_options, range_states)
    plot(P0_options, range_states[new_P = 1], label="Lower range", ylabel="New amount of Phosphorus",
         xlabel="Past amount of Phosphorus", left_margin = 5Plots.mm, size=(_cm(12), _cm(8)))
    plot!(P0_options, range_states[new_P = 2], label="Upper range")
    threshold = PathwayDiversity.get_root(1.3, 0.1)
    vline!([threshold], label="Tipping point", color="black", ls=:dash)
    hline!([threshold], label=false, color="black", ls=:dash)
    savefig("../output/range_states.pdf")
end

function _plot_number_option(n_possible_influx, P_options)
    plot(collect(P_options), n_possible_influx, label=false, ylabel="Number of options",
         xlabel="Amount of Phosphorus (x)", size=(_cm(12), _cm(8)))
    savefig("../output/number_options.png")

end

function _plot_bifurcation(roots, influx_options_root, lower_fold, upper_fold)
    # Each branch exists only on part of the influx range, bounded by the folds:
    # upper-stable above lower_fold, unstable between the folds, lower-stable below
    # upper_fold. Mask by influx value so it is independent of the sampling step.
    upper_branch  = influx_options_root .>= lower_fold
    middle_branch = lower_fold .<= influx_options_root .<= upper_fold
    lower_branch  = influx_options_root .<= upper_fold
    plot(influx_options_root[upper_branch], roots[3, upper_branch], lw=1.5)
    plot!(influx_options_root[middle_branch], roots[2, middle_branch], lw=1.5)
    plot!(influx_options_root[lower_branch], roots[1, lower_branch], lw=1.5,
          bottom_margin = 5Plots.mm, left_margin = 10Plots.mm,
          legend=false, ylabel = "Amount of phosphorus (x)", xlabel = "Influx (I)", size=(_cm(15), _cm(10)))
    savefig("../output/bifurcation.pdf")
end

function _clean_vector(vector)::Vector{Float64}
    first_string = findfirst(x -> x == "", vector)
    if first_string == nothing
        return Vector{Float64}(vector)
    end

    new_vector = Vector{Float64}(vector[begin:first_string-1])
    return new_vector
end
