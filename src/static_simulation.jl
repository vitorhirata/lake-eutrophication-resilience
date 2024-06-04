# run entropy for different initial conditions and time horizon. Used in distance to basin threshold analysis
function run_entropy(
    P_init_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizons::Vector{Float64},
)::String
    s = NamedDimsArray{(:P0, :time_horizon)}(zeros(length(P_init_options), length(time_horizons)))

    for (idx_P0, P0) in enumerate(P_init_options), (idx_time_horizon, time_horizon) in enumerate(time_horizons)
        number_decision::Int64 = floor(time_horizon / decision_step)
        s[idx_P0, idx_time_horizon] = _entropy(P0, influx, decision_step, number_decision)
        time_horizon == time_horizons[end] && P0 % 0.5 == 0 && println("Finished model for P0=$(P0)")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_distance_basin_"
    writedlm("$(base_filename)s.csv",  s, ',')
    return timestamp
end

# run entropy for different decision step and time horizon. Used in scalling analysis
function run_entropy(
    P_init_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizon::Float64,
    number_options::Vector{Int64},
)::String
    s = NamedDimsArray{(:number_options, :P0)}(zeros(length(number_options), length(P_init_options)))
    number_decision::Int64 = floor(time_horizon / decision_step)

    for (idx_P0, P0) in enumerate(P_init_options), (idx_number_option, number_option) in enumerate(number_options)
        s[idx_number_option, idx_P0] = _entropy(P0, influx, decision_step, number_decision, true, number_option)
        s[idx_number_option, idx_P0] /= (number_decision * log(number_option))
        number_option == number_options[end] && P0 % 1.0 == 0 && println("Finished model for P0=$(P0)")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_scale_initial_"
    writedlm("$(base_filename)s.csv",  s, ',')
    return timestamp
end

# run entropy for different initial conditions and maximum number of options. Used in decision scale analysis
function run_entropy(
    P0::Float64,
    influx::Float64,
    decision_steps::Vector{Float64},
    time_horizons::Vector{Float64},
)::String
    s = NamedDimsArray{(:decision_step, :time_horizon)}(zeros(length(decision_steps), length(time_horizons)))

    for (idx_step, step) in enumerate(decision_steps), (idx_time_horizon, time_horizon) in enumerate(time_horizons)
        number_decision::Int64 = floor(time_horizon / step)
        s[idx_step, idx_time_horizon] = _entropy(P0, influx, step, number_decision)
        time_horizon == time_horizons[end] && println("Finished model for P0=$(P0)")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_decision_scale_"
    writedlm("$(base_filename)s.csv",  s, ',')
    return timestamp
end

# evolve one time step and compute the number of options available for each P0 and influx.
function run_number_options_simulation(
    P0_options::Vector{Float64},
    influx_options::Vector{Float64},
    max_number_options::Int64 = 20,
    decision_step::Float64 = 5.0,
    deterministic::Bool = true
)
    result = NamedDimsArray{(:P0, :influx)}(zeros(length(P0_options), length(influx_options)))

    for (idx_P0, P0) in enumerate(P0_options), (idx_influx, influx) in enumerate(influx_options)
        new_P = _evolve_step(P0, influx, decision_step, deterministic)
        result[P0=idx_P0, influx=idx_influx] = number_possible_influx(new_P, max_number_options)
    end
    return result
end
