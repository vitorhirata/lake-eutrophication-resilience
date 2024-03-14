# run entropy for different initial conditions and time horizon. Used in distance to basin threshold analysis
function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizons::Vector{Float64},
)::NamedDimsArray
    s = NamedDimsArray{(:P0, :time_horizon)}(zeros(length(P0_options), length(time_horizons)))

    for (idx_P0, P0) in enumerate(P0_options), (idx_time_horizon, time_horizon) in enumerate(time_horizons)
        number_decision::Int64 = floor(time_horizon / decision_step)
        s[idx_P0, idx_time_horizon] = _entropy(P0, influx, decision_step, number_decision)
        time_horizon == time_horizons[end] && P0 % 0.5 == 0 && println("Finished model for P0=$(P0)")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_distance_basin_"
    writedlm("$(base_filename)s.csv",  s, ',')
    return s
end

# run entropy for different initial conditions and maximum number of options. Used in decision scale analysis
function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizon::Float64,
    number_options::Vector{Int64},
)::NamedDimsArray
    s = NamedDimsArray{(:number_options, :P0)}(zeros(length(number_options), length(P0_options)))
    number_decision::Int64 = floor(time_horizon / decision_step)

    for (idx_P0, P0) in enumerate(P0_options), (idx_number_option, number_option) in enumerate(number_options)
        s[idx_number_option, idx_P0] = _entropy(P0, influx, decision_step, number_decision, true, number_option)
        s[idx_number_option, idx_P0] /= (number_decision * log(number_option))
        number_option == number_options[end] && P0 % 1.0 == 0 && println("Finished model for P0=$(P0)")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_scale_initial_"
    writedlm("$(base_filename)s.csv",  s, ',')
    return s
end

# run entropy for different decision step and time horizon. Used in scalling analysis
function run_entropy(
    P0::Float64,
    influx::Float64,
    decision_steps::Vector{Float64},
    time_horizons::Vector{Float64},
)::NamedDimsArray
    s = NamedDimsArray{(:decision_step, :time_horizon)}(zeros(length(decision_steps), length(time_horizons)))

    for (idx_step, step) in enumerate(decision_steps), (idx_time_horizon, time_horizon) in enumerate(time_horizons)
        number_decision::Int64 = floor(time_horizon / step)
        s[idx_step, idx_time_horizon] = _entropy(P0, influx, step, number_decision)
        time_horizon == time_horizons[end] && println("Finished model for P0=$(P0)")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_decision_scale"
    writedlm("$(base_filename)s.csv",  s, ',')
    return s
end
