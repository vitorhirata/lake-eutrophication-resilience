function _plot_early_warning_signals(p, s, residuals, variance_ts, autocorr_ts, time_horizons, times,
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
    savefig("../output/early_warning_signal.png")
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
