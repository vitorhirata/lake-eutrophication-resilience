using PathwayDiversity

function scaling()
    influx = 0.1
    decision_step = 5.0
    time_horizon = 20.0
    P0_options = collect(0:0.5:3)
    number_options = collect(10:20:80)

    timestamp = PathwayDiversity.run_entropy(P0_options, influx, decision_step, time_horizon, number_options)
    PathwayDiversity.scaling(P0_options, number_options, timestamp)
end

function decision_scales()
    P0 = 0.1
    influx = 0.1
    decision_steps = [4.0, 6.0, 8.0]
    time_horizons = [8.0, 18.0, 24.0]

    timestamp = PathwayDiversity.run_entropy(P0, influx, decision_steps, time_horizons)
    PathwayDiversity.decision_scales(time_horizons, decision_steps, timestamp)
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
