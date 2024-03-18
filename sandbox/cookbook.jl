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
    number_options = 10

    n_possible_influx = PathwayDiversity.number_possible_influx.(P, number_options)
    plot(collect(P), n_possible_influx, label=false, ylabel="Number of options", xlabel="Amount of phosphorus (x)",
         guidefontsize=12)
    savefig("../output/number_options.png")
end

function scaling()
    influx = 0.1
    decision_step = 5.0
    time_horizon = 20.0

    P_init_options = collect(0:0.5:3)
    number_options = collect(10:20:80)

    s = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizon, number_options)
    _plot_scaling(s, P_init_options, number_options)
    return s
end

function decision_scales()
    P_init = 0.1
    influx = 0.1
    decision_steps = [4.0, 6.0, 8.0]
    time_horizons = [8.0, 18.0, 24.0]

    s = PathwayDiversity.run_entropy(P_init, influx, decision_steps, time_horizons)
    _plot_decision_scales(s, time_horizons, decision_steps)
    return s
end

function distance_basin_threshold()
    influx = 0.1
    decision_step = 5.0
    P_init_options = collect(0:0.05:3)
    time_horizons = [10.0, 20.0, 35.0]

    s = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizons)
    s = PathwayDiversity.normalize_pd(s)
    s_diff = PathwayDiversity.finite_difference(s, 0.05)

    threshold = PathwayDiversity.get_root(1.3, influx)
    distance_threshold = threshold .- P_init_options

    _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons)
    return s, distance_threshold, s_diff
end

function early_warning_signals()
    parameters = Dict(
        :P_init => 0.27,
        :times => 1:0.125:75,
        :time_horizons => [5.0, 10.0, 20.0],
        :decision_step => 5.0,
        :influx => 0.03,
        :influx_tax => 0.001,
    )
    timestamp = PathwayDiversity.run_scenario(; parameters...)
    PathwayDiversity.early_warning_signals(timestamp; parameters...)
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

function _plot_scaling(s, P_init_options, number_options)
    label = map(P_init -> "Initial condition = $(P_init)", P_init_options)
    label = reshape(label, (1,length(P_init_options)))

    plot(number_options, s, label=label, left_margin = 5Plots.mm, legend=:outerbottomright,
         size=(952,560), ylabel = "Pathway diversity", xlabel = "Maximum number of options")
    savefig("../output/scaling.png")
end

function _plot_decision_scales(s, time_horizons, decision_steps)
    label = map(decision_step -> "Time horizon = $(decision_step)", time_horizons)
    label = reshape(label, (1,length(time_horizons)))

    plot(decision_steps, s, label = label, legend=:topright, size=(952,560),
         ylabel = "Pathway diversity", xlabel = "Decision step", guidefontsize=12, left_margin = 10Plots.mm)
    savefig("../output/decision_scales.png")
end

function _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons)
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
          ylabel="Pathway diversity derivate", xlabel="Distance to threshold", left_margin = 10Plots.mm)
    vline!([0.0], label=false, color="black")

    plot(plt1, plt2, layout=(2,1), size=(1000,1200), guidefontsize=12)
    savefig("../output/distance_threshold.png")
end


