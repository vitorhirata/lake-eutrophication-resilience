using PathwayDiversity

function distance_basin_threshold()
    influx = 0.1
    decision_step = 5.0
    P0_options = collect(0:0.05:2.85)
    time_horizons = [7] * decision_step # [1, 2, 4, 7]

    timestamp = PathwayDiversity.run_entropy(P0_options, influx, decision_step, time_horizons)
    PathwayDiversity.distance_basin_threshold(P0_options, influx, time_horizons, decision_step, timestamp)
end

function states_distribution()
    P0_options = [0.5, 1.0, 1.5, 2.0]
    time_horizons = 35.0
    decision_steps = 5.0

    timestamp = PathwayDiversity.run_state_distribution(P0_options, time_horizons, decision_steps)
    PathwayDiversity.states_distribution(P0_options, time_horizons, decision_steps, timestamp)
end

function early_warning_signals()
    parameters = Dict(
        :P0 => 0.27,
        :times => 1:0.125:80,
        :time_horizons => [5.0, 10.0, 20.0],
        :decision_step => 5.0,
        :influx => 0.03,
        :influx_tax => 0.001,
    )
    timestamp = PathwayDiversity.run_time_series(; parameters...)
    PathwayDiversity.early_warning_signals(timestamp; parameters...)
end

function sensitivity()
    P0_options = collect(0:0.05:2.85)
    time_horizons = 30.0
    influx = 0.12
    scenarios = [
        Dict(:name => "Default", :decision_step => 5.0),
        Dict(:name => "1. Reduced decisions", :decision_step => 15.0),
        Dict(:name => "2. Allow restrictive option", :decision_step => 5.0, :minimum_influx => 0.02, :maximum_influx => 0.32),
        Dict(:name => "3. Probability based on change", :decision_step => 5.0, :method => "state_change"),
        Dict(:name => "4. Minimal change", :decision_step => 5.0, :method => "closer_more_likely"),
        Dict(:name => "5. Favour lower influx", :decision_step => 5.0, :method => "favor_positive"),
        #Dict(:name => "Unstable decisions", :decision_step => 5.0, :method => "further_more_likely"),
        #Dict(:name => "Favor negative", :decision_step => 5.0, :method => "favor_negative"),
    ]

    timestamp = PathwayDiversity.run_sensitivity(P0_options, influx, time_horizons, scenarios)
    PathwayDiversity.sensitivity(P0_options, influx, time_horizons, scenarios, timestamp, true)
end
