function _plot_early_warning_signals(timestamp, s, variance_ts, autocorr_ts, distance_thresholds,
    time_horizons, times, tipping_points)

    label = reshape(map(time_horizon -> "Time horizon = $(time_horizon)", time_horizons), (1, size(s, :time_horizon)))
    xticks = 0:25:length(times)
    xlims = (0, times[end]+1)

    plt1 = plot(collect(times), distance_thresholds, label=false, ylabel="Distance to threshold")
    vline!([times[tipping_points[:p]]], label="Threshold", color="black", lw=2, xticks=xticks, xlims=xlims)

    plt2 = plot(collect(times), s, label=label, ylabel="Pathway diversity", legend=:bottomleft)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims)
    for horizon in 1:size(s, :time_horizon)
        scatter!([times[tipping_points[:s][horizon]]], [s[time=tipping_points[:s][horizon], time_horizon=horizon]],
                 label=false, markerstrokewidth=0, color=horizon)
    end

    new_time = times[1 + length(times) - length(variance_ts):end]
    plt3 = plot(collect(new_time), variance_ts, label=false, ylabel="Variance")
    scatter!([new_time[tipping_points[:var]]], [variance_ts[tipping_points[:var]]],
             label="Kendall-τ > 0.56", markerstrokewidth=0, color=1)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims)

    new_time = times[1 + length(times) - length(autocorr_ts):end]
    plt4 = plot(collect(new_time), autocorr_ts, label=false, ylabel="Autocorrelation")
    scatter!([new_time[tipping_points[:autocorr]]], [autocorr_ts[tipping_points[:autocorr]]],
             label=false, markerstrokewidth=0, color=1, xlabel="Time (year)")
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2, xticks=xticks, xlims=xlims)

    plot(plt1, plt2, plt3, plt4, layout=(4,1), size=(1000,1000), guidefontsize=12, left_margin = 10Plots.mm)
    savefig("../output/$(timestamp)_early_warning_signal.png")
end

function _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons, timestamp, peaks_idx, one_plot)
    label = reshape(map(t_horizon -> "Time horizon = $(t_horizon)", time_horizons), (1,size(s, :time_horizon)))
    peak_values = [s_diff[P0=peak_idx, time_horizon=idx] for (idx,peak_idx) in enumerate(peaks_idx)]
    xlims = (distance_threshold[end] - 0.1, distance_threshold[begin] + 0.1)

    if one_plot
        plt1 = plot(distance_threshold, s, label=label, legend=:left)
        vline!([0.0], label=false, color="black", xflip = true, ylabel="Pathway diversity", xlims=xlims)

        plt2 = plot(distance_threshold[2:end], s_diff, label=label, legend=:topleft, xlabel="Distance to threshold")
        scatter!(distance_threshold[2:end][peaks_idx], peak_values, label=false, markerstrokewidth=0, color=[1,2,3,4])
        vline!([0.0], label=false, color="black", xflip = true, ylabel="Pathway diversity derivative", xlims=xlims)

        plot(plt1, plt2, layout=(2,1), size=(1000,1200), guidefontsize=12, left_margin = 10Plots.mm)
        savefig("../output/$(timestamp)_distance_threshold.png")
    else
        plot(distance_threshold, s, label=label, legend=:left, xlabel="Distance to threshold")
        vline!([0.0], label=false, color="black", xflip = true, ylabel="Pathway diversity",
               guidefontsize=12, left_margin = 10Plots.mm, size=(1000,600), xlims=xlims)
        savefig("../output/$(timestamp)_distance_threshold.png")

        plot(distance_threshold[2:end], s_diff, label=label, legend=:topleft, xlabel="Distance to threshold")
        scatter!(distance_threshold[2:end][peaks_idx], peak_values, label=false, markerstrokewidth=0, color=[1,2,3,4])
        vline!([0.0], label=false, color="black", xflip = true, ylabel="Pathway diversity derivative",
               guidefontsize=12, left_margin = 10Plots.mm, size=(1000,600), xlims=xlims)
        savefig("../output/$(timestamp)_distance_threshold_derivative.png")
    end
end

function _plot_sensitivity(s, distance_threshold, timestamp, scenarios, relative = false)
    if relative
        first_idx = 2
        plot(distance_threshold, s[type=first_idx], label="$(scenarios[first_idx][:name])", legend=:topright,
             ylabel="Pathway diversity", left_margin = 10Plots.mm, size=(1000,600),
             alpha=0.9, guidefontsize=12, xflip = true, xlabel="Distance to threshold", ylims=[0.67, 1.3])
        hline!([1.0], label=false, color="black")
    else
        first_idx = 1
        plot(distance_threshold, s[type=first_idx], label="$(scenarios[first_idx][:name])", legend=:topright,
             ylabel="Pathway diversity", left_margin = 10Plots.mm, size=(1000,600), color="black", lw=1.5,
             alpha=0.9, guidefontsize=12, xflip = true, xlabel="Distance to threshold")
    end

    vline!([0.0], label=false, color="black")
    for (idx_scenario, scenario) in enumerate(scenarios[first_idx+1:end])
        plot!(distance_threshold, s[type=idx_scenario+first_idx], label="$(scenario[:name])",
              color=idx_scenario+first_idx-1, alpha=0.9)
    end
    savefig("../output/$(timestamp)_sensitivity.png")
end

function _plot_states_distribution(P0_options, n_decision, timestamp)
    plot_array = Plots.Plot[]
    labels = reshape(map(decision -> "Decision $(decision)", 1:n_decision), 1, n_decision)

    for (idx, P0) in enumerate(P0_options)
        states = readdlm("../output/$(timestamp)_state_distribution_$(idx).csv", ',')
        states = [_clean_vector(row) for row in eachrow(states)]
        plt1 = violin(labels, states[2:end], label=false, title="Distance to threshold: $(round(P0, digits=2))",
                      ylim=(0.0, 2.5), ylabel="State", color=idx)
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
