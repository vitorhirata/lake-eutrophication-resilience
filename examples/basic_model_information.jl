using PathwayDiversity

# Plot the bifurcation diagram over influx values from 0.00 to 0.30.
# Shows stable and unstable equilibria as a function of influx.
function bifurcation()
    influx_options_root = collect(range(0.00, 0.3,  step = 0.001))

    PathwayDiversity.bifurcation(influx_options_root)
end

# Plot available influx options as a function of phosphorus level.
function number_options()
    P_options = 0.0:0.05:4.0
    max_number_options = 10

    PathwayDiversity.number_options(P_options, max_number_options)
end

# Simulate and plot influx option counts for a (P0, influx) grid.
# Uses up to 20 options and covers P0 from 0 to 3.5.
function number_options_simulation()
    P0_options = collect(0:0.01:3.5)
    influx_options = collect(range(0.00, 0.3, step = 0.001))
    max_number_options = 20

    PathwayDiversity.number_options_simulation(P0_options, influx_options, max_number_options)
end

# Simulate and plot min/max reachable states after one decision step.
# Covers P0 from 0 to 3.5 with up to 10 available influx options.
function range_state()
    P0_options = collect(0:0.01:3.5)
    max_number_options = 10

    PathwayDiversity.range_states_simulation(P0_options, max_number_options)
end
