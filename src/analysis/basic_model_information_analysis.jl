function bifurcation(influx_options_root::Vector{Float64})
    roots = zeros(3, length(influx_options_root))

    for (index, influx) in enumerate(influx_options_root)
        roots[1, index] = PathwayDiversity.get_root(0.1, influx)
        roots[2, index] = PathwayDiversity.get_root(1.3, influx)
        roots[3, index] = PathwayDiversity.get_root(2.7, influx)
    end
    _plot_bifurcation(roots, influx_options_root)
end

function number_options(P_options::StepRangeLen{Float64}, max_number_options::Int64)
    n_possible_influx = PathwayDiversity.number_possible_influx.(P_options, max_number_options)
    _plot_number_option(n_possible_influx, P_options)
end

function number_options_simulation(
    P0_options::Vector{Float64}, influx_options::Vector{Float64}, max_number_options::Int64
)
    n_options = PathwayDiversity.run_number_options_simulation(P0_options, influx_options, max_number_options)
    _plot_number_options_simulation(P0_options, influx_options, n_options)
end

function range_states_simulation(P0_options::Vector{Float64}, max_number_options::Int64)
    range_states = PathwayDiversity.run_range_states_simulation(P0_options, max_number_options)
    _plot_range_states(P0_options, range_states)
end
