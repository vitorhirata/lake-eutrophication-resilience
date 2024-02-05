using Revise
using Infiltrator
using PathwayDiversity
using Plots
using NamedDims

function all()
    bifurcation()
    number_options()
    scaling()
    decision_scales()
    distance_basin_threshold()
    early_warning_signals()
end

function bifurcation()
    influx_options_root = influx_options_root2 = range(0.00, 0.3,  step = 0.001) |> collect
    roots = zeros(3, length(influx_options_root)  )

    for (index, influx) in enumerate(influx_options_root)
        roots[1, index] = PathwayDiversity.get_root(0.1, influx)
        roots[2, index] = PathwayDiversity.get_root(1.3, influx)
        roots[3, index] = PathwayDiversity.get_root(2.7, influx)
    end
    _plot_bifurcation(roots, influx_options_root)
    return roots
end

function number_options()
    P = 0.0:0.05:4.0
    number_options = 10:5:30

    plot(ylabel="Number of options", xlabel="Amount of phosphorus (P)")
    for number_option in number_options
        n_possible_influx = PathwayDiversity.number_possible_influx.(P, number_option)
        plot!(collect(P), n_possible_influx, label="Maximum number of options = $(number_option)")
    end
    savefig("../output/number_options.png")
end

function scaling()
    influx = 0.1
    decision_step = 5.0
    time_horizon = 20.0

    P_init_options = collect(0:0.5:3)
    number_options = collect(10:20:80)

    s_final = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizon, number_options)
    _plot_scaling(s_final, P_init_options, number_options)
    return s_final
end

function decision_scales()
    P_init = 0.1
    influx = 0.1
    decision_steps = [4.0, 6.0, 8.0]
    time_horizons = [8.0, 18.0, 24.0]

    s_final = PathwayDiversity.run_entropy(P_init, influx, decision_steps, time_horizons)
    _plot_decision_scales(s_final, time_horizons, decision_steps)
    return s_final
end

function distance_basin_threshold()
    influx = 0.1
    decision_step = 5.0
    P_init_options = collect(0:0.1:2.4)
    time_horizons = collect(decision_step:decision_step:35.0)

    s_final = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizons)

    _plot_time_inital_state(s_final, P_init_options, time_horizons)

    threshold = PathwayDiversity.get_root(1.3, influx)
    _plot_distance_threshold(s_final, P_init_options, time_horizons, threshold)
    return s_final, threshold
end

function early_warning_signals()
    P_init = 0.02
    t_max = 150
    step = 0.125
    time_horizon = 20.0
    decision_step = 5.0
    influx = 0.03
    influx_taxes = [0.0, 0.0001, 0.0005, 0.001, 0.0015, 0.002]

    result = PathwayDiversity.run_scenarios(P_init, influx, influx_taxes, t_max, step, decision_step, time_horizon)

    times = collect(1:step:t_max)
    residuals = PathwayDiversity.detrend(result[type=1], times)

    # Compute variance
    variance_time_step = 5
    variance_idx_step::Int64 = variance_time_step ÷ step
    variance_ts = PathwayDiversity.compute_variance(residuals, variance_idx_step)

    # Compute autocorrelation
    autocorr_time_step = 25
    autocorr_idx_step::Int64 = autocorr_time_step ÷ step
    autocorr_ts = PathwayDiversity.compute_autocorrelation(residuals, autocorr_idx_step)

    _plot_early_warning_signals(result, residuals, variance_ts, autocorr_ts, influx_taxes,
                                step, t_max, variance_time_step, autocorr_time_step)

    return result, residuals, variance_ts, autocorr_ts
end

function _plot_bifurcation(roots, influx_options_root)
    lim1 = 19
    lim2 = 174
    plot(influx_options_root[lim1:end], roots[3, lim1:end], label="Eutrophicated stable state")
    plot!(influx_options_root[lim1:lim2], roots[2, lim1:lim2], label="Instable state")
    plot!(influx_options_root[1:lim2], roots[1, 1:lim2], label="Clean stable state")

    hline!([0.4], label="Number of options drop", legend=:bottomright, size=(952,560),
           ylabel = "Fixed points (P*)", xlabel = "Influx (I)")
    savefig("../output/bifurcation.png")
end

function _plot_scaling(s_final, P_init_options, number_options)
    label = map(P_init -> "Initial condition = $(P_init)", P_init_options)
    label = reshape(label, (1,length(P_init_options)))

    plot(number_options, transpose(s_final), label=label, left_margin = 5Plots.mm, legend=:outerbottomright,
         size=(952,560), ylabel = "Pathway diversity", xlabel = "Maximum number of options")
    savefig("../output/scaling.png")
end

function _plot_decision_scales(s_final, time_horizons, decision_steps)
    label = map(decision_step -> "Time horizon = $(decision_step)", time_horizons)
    label = reshape(label, (1,length(time_horizons)))

    plot(decision_steps, s_final, label = label, legend=:topright, size=(952,560),
         ylabel = "Pathway diversity", xlabel = "Decision step")
    savefig("../output/decision_scales.png")
end

function _plot_time_inital_state(s_final, P_init_options, time_horizons)
    label = map(P_init -> "Initial state (P0) = $(P_init)", P_init_options)
    selected_index = [1, 4, 7, 10, 13, 16, 19, 25]

    s_final_filtered = stack([s_final[i, :] for i in selected_index], dims=1)
    label = [label[i] for i in selected_index]
    label = reshape(label, (1,length(selected_index)))

    plot(time_horizons, transpose(s_final_filtered), label=label, legend=:topleft, size=(952,560),
         ylabel = "Pathway diversity", xlabel = "Time horizon")
    savefig("../output/time_initial_state.png")
end

function _plot_distance_threshold(s_final, P_init_options, time_horizons, threshold)
    label = map(t_horizon -> "Time horizon = $(t_horizon)", time_horizons)
    selected_index = [1, 2, 3, 5, 7]

    s_final_filtered = stack([s_final[:, i] for i in selected_index], dims=1)
    label = [label[i] for i in selected_index]
    label = reshape(label, (1,length(selected_index)))
    distance_threshold = threshold .- P_init_options

    plot(distance_threshold, transpose(s_final_filtered), label=label, legend=:topleft, size=(952,560),
         ylabel="Pathway diversity", xlabel="Distance to threshold")
    vline!([0.0], label=false)

    savefig("../output/distance_threshold.png")
end

function _plot_early_warning_signals(result, residuals, variance_ts, autocorr_ts, influx_taxes, step, t_max,
        variance_time_step, autocorr_time_step
)
    label = map(influx_tax -> "Influx_tax = $(influx_tax)", influx_taxes)
    xticks = 0:25:length(1:step:t_max)
    xlims = (0, t_max+1)

    selected_index = [1, 2, 3, 4, 5, 6]
    label = reshape(label[selected_index], (1,length(selected_index)))
    p = result[type=1, influx_tax=selected_index]
    s = result[type=2, influx_tax=selected_index]
    residual = residuals[influx_tax=selected_index]

    variance_idx_step::Int64 = variance_time_step ÷ step
    autocorr_idx_step::Int64 = autocorr_time_step ÷ step
    variance = variance_ts[time=(variance_idx_step+1):end, influx_tax=selected_index]
    autocorr = autocorr_ts[time=(autocorr_idx_step+1):end, influx_tax=selected_index]

    plt1 = plot(collect(1:step:t_max), p, label=label, xticks=xticks, ylabel="Amount of Phosphorus",
                xlims=xlims, left_margin = 5Plots.mm)
    plt2 = plot(collect(1:step:t_max), residual, label=false, ylabel="Residuals",
                xticks=xticks, xlims=xlims, left_margin = 5Plots.mm)
    plt3 = plot(collect(1:step:t_max), s, label=false, ylabel="Pathway diversity",
                xticks=xticks, xlims=xlims, left_margin = 5Plots.mm)
    plt4 = plot(collect((variance_time_step+1):step:t_max), variance, label=false, xticks=xticks,
                ylabel="Variance", xlims=xlims, left_margin = 10Plots.mm)
    plt5 = plot(collect((autocorr_time_step+1):step:t_max), autocorr, label=false, xticks=xticks,
                ylabel="Aucorrelation", xlabel="Time (year)", xlims=xlims, left_margin = 10Plots.mm)

    plot(plt1, plt2, plt3, plt4, plt5, layout=(5,1), legend=:topleft, size=(1250,900), guidefontsize=10)
    savefig("../output/early_warning_signal.png")
end

