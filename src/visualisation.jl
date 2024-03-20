function _plot_early_warning_signals(timestamp, p, s, residuals, variance_ts, autocorr_ts, time_horizons, times,
        variance_time_step, autocorr_time_step, tipping_points, kendall_tau, include_residual = false
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
    plt3 = plot(collect(times), s, label=label, ylabel="Pathway diversity",
                xticks=xticks, xlims=xlims, left_margin = 5Plots.mm, legend=:bottomleft)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)
    for horizon in selected_index
        scatter!([times[tipping_points[:s][horizon]]], [s[time=tipping_points[:s][horizon], time_horizon=horizon]],
                 label="Kendall-τ=$(kendall_tau[:s][horizon])", markerstrokewidth=0, color=horizon)
    end
    plt4 = plot(collect((variance_time_step+1):step(times):times[end]), variance, label=false, xticks=xticks,
                ylabel="Variance", xlims=xlims, left_margin = 10Plots.mm)
    scatter!([times[tipping_points[:var]]], [variance[tipping_points[:var]-(variance_idx_step+1)]],
             label="Kendall-τ=$(kendall_tau[:var])", markerstrokewidth=0, color=1)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)
    plt5 = plot(collect((autocorr_time_step+1):step(times):times[end]), autocorr, label=false, xticks=xticks,
                ylabel="Autocorrelation", xlabel="Time (year)", xlims=xlims, left_margin = 10Plots.mm)
    scatter!([times[tipping_points[:autocorr]]], [autocorr[tipping_points[:autocorr]-(autocorr_idx_step+1)]],
             label="Kendall-τ=$(kendall_tau[:autocorr])", markerstrokewidth=0, color=1)
    vline!([times[tipping_points[:p]]], label=false, color="black", lw=2)

    if include_residual
        plot(plt1, plt2, plt3, plt4, plt5, layout=(5,1), size=(1000,900), guidefontsize=12)
    else
        plot(plt1, plt3, plt4, plt5, layout=(4,1), size=(1000,900), guidefontsize=12)
    end
    savefig("../output/$(timestamp)_early_warning_signal.png")
end

function _plot_kendall_taus(kendall_taus, time_horizons)
    labels = ["variance", "autocorrelation"]
    for time_horizon in time_horizons
        push!(labels, "PD_$(time_horizon)")
    end
    xticks = (1:(2+length(time_horizons)), labels)

    plt1 = boxplot(parent(kendall_taus[type=1]), label=false, ylabel="Tipping point time")
    plot!(xticks=xticks)

    plt2 = boxplot(parent(kendall_taus[type=2]), label=false, ylabel="Kendall-τ value")
    plot!(xticks=xticks)

    plot(plt1, plt2, layout=(2,1), size=(600,500), guidefontsize=12)
    savefig("../output/kendall_tau.png")
end

function _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons, timestamp)
    label = map(t_horizon -> "Time horizon = $(t_horizon)", time_horizons)
    selected_index = [1, 2, 3]

    s = s[time_horizon=selected_index]
    s_diff = s_diff[time_horizon=selected_index]
    label = [label[i] for i in selected_index]
    label = reshape(label, (1,length(selected_index)))

    plt1 = plot(distance_threshold, s, label=label, legend=:left, xflip = true,
          ylabel="Pathway diversity", left_margin = 10Plots.mm)
    vline!([0.0], label=false, color="black")
    plt2 = plot(distance_threshold[2:end], s_diff, label=label, legend=:left, color=[1 2 3], xflip = true,
          ylabel="Pathway diversity derivative", xlabel="Distance to threshold", left_margin = 10Plots.mm)
    vline!([0.0], label=false, color="black")

    plot(plt1, plt2, layout=(2,1), size=(1000,1200), guidefontsize=12)
    savefig("../output/$(timestamp)_distance_threshold.png")
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
