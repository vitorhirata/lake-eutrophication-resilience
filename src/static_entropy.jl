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
    end
    return s
end

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
    end
    return s
end

function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizon::Float64,
    number_options::Vector{Int64}
)::NamedDimsArray
    s = NamedDimsArray{(:number_options, :P0)}(zeros(length(number_options), length(P0_options)))
    number_decision::Int64 = floor(time_horizon / decision_step)

    for (idx_P0, P0) in enumerate(P0_options), (idx_number_option, number_option) in enumerate(number_options)
        s[idx_number_option, idx_P0] = _entropy(P0, influx, decision_step, number_decision, true, number_option)
        s[idx_number_option, idx_P0] /= (number_decision * log(number_option))
    end
    return s
end
