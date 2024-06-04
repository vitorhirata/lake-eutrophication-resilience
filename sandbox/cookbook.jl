using Revise
using Infiltrator
using PathwayDiversity

function all()
    bifurcation()
    number_options()
    number_options_simulation()
    range_state()
    scaling()
    decision_scales()
    distance_basin_threshold()
    early_warning_signals()
end

function bifurcation()
    influx_options_root = collect(range(0.00, 0.3,  step = 0.001))

    PathwayDiversity.bifurcation(influx_options_root)
end

function number_options()
    P_options = 0.0:0.05:4.0
    max_number_options = 10

    PathwayDiversity.number_options(P_options, max_number_options)
end

function number_options_simulation()
    P0_options = collect(0:0.01:3.5)
    influx_options = collect(range(0.00, 0.3, step = 0.001))
    max_number_options = 20

    PathwayDiversity.number_options_simulation(P0_options, influx_options, max_number_options)
end

function range_state()
    P0_options = collect(0:0.01:3.5)
    max_number_options = 10

    PathwayDiversity.range_states_simulation(P0_options, max_number_options)
end

function scaling()
    influx = 0.1
    decision_step = 5.0
    time_horizon = 20.0
    P_init_options = collect(0:0.5:3)
    number_options = collect(10:20:80)

    timestamp = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizon, number_options)
    PathwayDiversity.scaling(P_init_options, number_options, timestamp)
end

function decision_scales()
    P_init = 0.1
    influx = 0.1
    decision_steps = [4.0, 6.0, 8.0]
    time_horizons = [8.0, 18.0, 24.0]

    timestamp = PathwayDiversity.run_entropy(P_init, influx, decision_steps, time_horizons)
    PathwayDiversity.decision_scales(time_horizons, decision_steps, timestamp)
end

function distance_basin_threshold()
    influx = 0.1
    decision_step = 5.0
    P_init_options = collect(0:0.05:3)
    time_horizons = [1, 2, 4, 7] * decision_step

    timestamp = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizons)
    PathwayDiversity.distance_basin_threshold(P_init_options, influx, time_horizons, timestamp)
end

function early_warning_signals()
    parameters = Dict(
        :P_init => 0.27,
        :times => 1:0.125:80,
        :time_horizons => [5.0, 10.0, 20.0],
        :decision_step => 5.0,
        :influx => 0.03,
        :influx_tax => 0.001,
    )
    timestamp = PathwayDiversity.run_scenario(; parameters...)
    PathwayDiversity.early_warning_signals(timestamp; parameters...)
end
