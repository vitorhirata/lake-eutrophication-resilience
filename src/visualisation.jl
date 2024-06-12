function _plot_early_warning_signals(timestamp, p, s, residuals, variance_ts, autocorr_ts, influx_ts, thresholds,
    time_horizons, times, variance_time_step, autocorr_time_step, tipping_points, kendall_tau, include_residual = false
)
    label = map(time_horizon -> "Time horizon = $(time_horizon)", time_horizons)
    xticks = 0:25:length(times)
    xlims = (0, times[end]+1)

    selected_index = [1, 2, 3]
    label = reshape(label[selected_index], (1,length(selected_index)))
    s = s[time_horizon=selected_index]

    variance_idx_step::Int64 = variance_time_step ÷ step(times)
    autocorr_idx_step::Int64 = autocorr_time_step ÷ step(times)
    variance = variance_ts[(variance_idx_step+1):end]
    autocorr = autocorr_ts[(autocorr_idx_step+1):end]

    plt1 = plot(collect(times), p, label=false, xticks=xticks, ylabel="Amount of Phosphorus",
                xlims=xlims, left_margin = 5Plots.mm)
    vline!([times[tipping_points[:p]]], label="Tipping point", color="black", lw=2)
    plt2 = plot(collect(times), residuals, label=false, ylabel="Residuals",
                xticks=xticks, xlims=xlims, left_margin = 5Plots.mm)
    plt3 = plot(collect(times), thresholds, label=false, ylabel="Distance to threshold",
                xticks=xticks, xlims=xlims, left_margin = 5Plots.mm)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)
    plt4 = plot(collect(times), s, label=label, ylabel="Pathway diversity",
                xticks=xticks, xlims=xlims, left_margin = 5Plots.mm, legend=:bottomleft)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)
    for horizon in selected_index
        scatter!([times[tipping_points[:s][horizon]]], [s[time=tipping_points[:s][horizon], time_horizon=horizon]],
                 label=false, markerstrokewidth=0, color=horizon)
    end
    plot5 = plot(collect((variance_time_step+1):step(times):times[end]), variance, label=false, xticks=xticks,
                ylabel="Variance", xlims=xlims, left_margin = 10Plots.mm)
    scatter!([times[tipping_points[:var]]], [variance[tipping_points[:var]-(variance_idx_step+1)]],
             label=false, markerstrokewidth=0, color=1)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)
    plt6 = plot(collect((autocorr_time_step+1):step(times):times[end]), autocorr, label=false, xticks=xticks,
                ylabel="Autocorrelation", xlabel="Time (year)", xlims=xlims, left_margin = 10Plots.mm)
    scatter!([times[tipping_points[:autocorr]]], [autocorr[tipping_points[:autocorr]-(autocorr_idx_step+1)]],
             label=false, markerstrokewidth=0, color=1)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)

    if include_residual
        plot(plt1, plt2, plt3, plt4, plot5, plt6, layout=(6,1), size=(1000,1000), guidefontsize=12)
    else
        plot(plt1, plt3, plt4, plot5, plt6, layout=(5,1), size=(1000,1000), guidefontsize=12)
    end
    savefig("../output/$(timestamp)_early_warning_signal.png")
end

function _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons, timestamp, peaks_idx)
    label = map(t_horizon -> "Time horizon = $(t_horizon)", time_horizons)
    selected_index = [1, 2, 3, 4]

    s = s[time_horizon=selected_index]
    s_diff = s_diff[time_horizon=selected_index]
    label = [label[i] for i in selected_index]
    label = reshape(label, (1,length(selected_index)))
    peak_values = [s_diff[P0=peak_idx, time_horizon=idx] for (idx,peak_idx) in enumerate(peaks_idx)]

    plot(distance_threshold, s, label=label, legend=:left, xflip = true,
         ylabel="Pathway diversity", left_margin = 10Plots.mm, size=(1000,600), guidefontsize=12)
    vline!([0.0], label=false, color="black", xlabel="Distance to threshold")
    savefig("../output/$(timestamp)_distance_threshold.png")

    plot(distance_threshold[2:end], s_diff, label=label, legend=:left, color=[1 2 3 4], xflip = true,
          ylabel="Pathway diversity derivative", xlabel="Distance to threshold", left_margin = 10Plots.mm)
    scatter!(distance_threshold[2:end][peaks_idx], peak_values, label=false, markerstrokewidth=0, color=[1, 2, 3, 4])
    vline!([0.0], label=false, color="black", size=(1000,600), guidefontsize=12)
    savefig("../output/$(timestamp)_distance_threshold_derivative.png")
end

function _plot_states_distribution(P0_options, n_decision, timestamp)
    plot_array = Plots.Plot[]
    labels = reshape(map(decision -> "Decision $(decision)", 1:n_decision), 1, n_decision)

    for (idx, P0) in enumerate(P0_options)
        states = readdlm("../output/$(timestamp)_state_distribution_$(idx).csv", ',')
        states = [_clean_vector(row) for row in eachrow(states)]
        plt1 = violin(labels, states[2:end], label=false, title="Initial State $(P0)", ylabel="State", color=idx,
                      ylim=(0.0, 2.5))
        push!(plot_array,plt1)
    end
    plot(plot_array..., layout=(length(P0_options), 1), size=(800,200*length(P0_options)), guidefontsize=12)
    savefig("../output/$(timestamp)_state_distribution.png")
end

function _plot_decision_scales(s, time_horizons, decision_steps, timestamp)
    label = map(decision_step -> "Time horizon = $(decision_step)", time_horizons)
    label = reshape(label, (1,length(time_horizons)))

    plot(decision_steps, s, label = label, legend=:topright, size=(952,560),
         ylabel = "Pathway diversity", xlabel = "Decision step", guidefontsize=12, left_margin = 10Plots.mm)
    savefig("../output/$(timestamp)_decision_scales.png")
end


function _plot_scaling(s, P_init_options, number_options, timestamp)
    label = map(P_init -> "Initial condition = $(P_init)", P_init_options)
    label = reshape(label, (1,length(P_init_options)))

    plot(number_options, s, label=label, left_margin = 5Plots.mm, legend=:outerbottomright,
         size=(952,560), ylabel = "Pathway diversity", xlabel = "Maximum number of options")
    savefig("../output/$(timestamp)_scaling.png")
end

function _plot_number_options_simulation(P0_options, influx_options, n_options)
    heatmap(influx_options, P0_options, n_options, label=false, ylabel="Amount of phosphorus", xlabel="Influx",
            guidefontsize=12, colorbar_title="Future number of options", left_margin = 5Plots.mm, size=(952,560))
    savefig("../output/number_options_simulation.png")
end

function _plot_range_states(P0_options, range_states)
    plot(P0_options, range_states[new_P = 1], label="Lower range", ylabel="New amount of phosphorus",
         xlabel="Past amount of phosphorus", guidefontsize=12, left_margin = 5Plots.mm, size=(952,560))
    plot!(P0_options, range_states[new_P = 2], label="Upper range")
    threshold = PathwayDiversity.get_root(1.3, 0.1)
    vline!([threshold], label="Tipping point", color="black", ls=:dash)
    hline!([threshold], label=false, color="black", ls=:dash)
    savefig("../output/range_states.png")
end

function _plot_number_option(n_possible_influx, P_options)
    plot(collect(P_options), n_possible_influx, label=false, ylabel="Number of options",
         xlabel="Amount of phosphorus (x)", guidefontsize=12)
    savefig("../output/number_options.png")

end

function _plot_bifurcation(roots, influx_options_root)
    lim1 = 19
    lim2 = 174
    plot(influx_options_root[lim1:end], roots[3, lim1:end], label="Eutrophicated stable state")
    plot!(influx_options_root[lim1:lim2], roots[2, lim1:lim2], label="Instable state")
    plot!(influx_options_root[1:lim2], roots[1, 1:lim2], label="Clean stable state")

    hline!([0.4], label="Number of options drop", legend=:bottomright, size=(952,560),
           ylabel = "Fixed points (P*)", xlabel = "Influx (I)", left_margin = 10Plots.mm)
    savefig("../output/bifurcation.png")
end

function _clean_vector(vector)::Vector{Float64}
    first_string = findfirst(x -> x == "", vector)
    if first_string == nothing
        return Vector{Float64}(vector)
    end

    new_vector = Vector{Float64}(vector[begin:first_string-1])
    return new_vector
end
