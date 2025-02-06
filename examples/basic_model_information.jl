using PathwayDiversity

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
