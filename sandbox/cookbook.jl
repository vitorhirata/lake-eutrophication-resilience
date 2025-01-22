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
    states_distribution()
    sensitivity()
    distance_basin_threshold()
    early_warning_signals()
    new_pathway_diversity_computation()
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

function states_distribution()
    P0_options = [0.663, 0.913, 1.413, 1.663]
    time_horizons = 35.0
    decision_steps = 5.0

    timestamp = PathwayDiversity.run_state_distribution(P0_options, time_horizons, decision_steps)
    PathwayDiversity.states_distribution(P0_options, time_horizons, decision_steps, timestamp)
end

function sensitivity()
    P0_options = collect(0:0.05:2.85)
    time_horizons = 30.0
    influx = 0.12
    scenarios = [
        Dict(:name => "Default", :decision_step => 5.0),
        Dict(:name => "1. Reduced decisions", :decision_step => 15.0),
        Dict(:name => "2. Restrictive option", :decision_step => 5.0, :minimum_influx => 0.02, :maximum_influx => 0.32),
        Dict(:name => "3. Probability based on change", :decision_step => 5.0, :method => "state_change"),
        Dict(:name => "4. Minimal change", :decision_step => 5.0, :method => "closer_more_likely"),
        Dict(:name => "5. Favor positive", :decision_step => 5.0, :method => "favor_positive"),
        #Dict(:name => "Unstable decisions", :decision_step => 5.0, :method => "further_more_likely"),
        #Dict(:name => "Favor negative", :decision_step => 5.0, :method => "favor_negative"),
    ]

    timestamp = PathwayDiversity.run_sensitivity(P0_options, influx, time_horizons, scenarios)
    PathwayDiversity.sensitivity(P0_options, influx, time_horizons, scenarios, timestamp, true)
end

function distance_basin_threshold()
    influx = 0.1
    decision_step = 5.0
    P_init_options = collect(0:0.05:2.85)
    time_horizons = [1, 2, 4, 7] * decision_step

    timestamp = PathwayDiversity.run_entropy(P_init_options, influx, decision_step, time_horizons)
    PathwayDiversity.distance_basin_threshold(P_init_options, influx, time_horizons, decision_step, timestamp)
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
    timestamp = PathwayDiversity.run_time_series(; parameters...)
    PathwayDiversity.early_warning_signals(timestamp; parameters...)
end

function new_pathway_diversity_computation()
    influx = 0.15
    decision_step = 5.0
    P0 = 0.5
    number_decision = 2

    old_entropy = PathwayDiversity.entropy(P0, influx, decision_step, number_decision)
    new_entropy = PathwayDiversity.new_entropy(P0, influx, decision_step, number_decision)

    print("Old method: $(old_entropy) \nNew method: $(new_entropy)")
end
