function run_entropy(
    P0::Float64,
    influx::Float64,
    decision_steps::Vector{Float64},
    time_horizons::Vector{Float64},
)::Matrix{Float64}
    s_final = zeros(length(decision_steps), length(time_horizons))

    for (idx_decision_step, decision_step) in enumerate(decision_steps)
        for (idx_time_horizon, time_horizon) in enumerate(time_horizons)
            number_decision::Int64 = floor(time_horizon / decision_step)
            s_final[idx_decision_step, idx_time_horizon] = _entropy(P0, influx, decision_step, number_decision)
        end
    end
    return s_final
end

function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizons::Vector{Float64},
)::Matrix{Float64}
    s_final = zeros(length(P0_options), length(time_horizons))

    for (idx_P0, P0) in enumerate(P0_options)
        for (idx_time_horizon, time_horizon) in enumerate(time_horizons)
            number_decision::Int64 = floor(time_horizon / decision_step)
            s_final[idx_P0, idx_time_horizon] = _entropy(P0, influx, decision_step, number_decision)
        end
    end
    return s_final
end

function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    time_horizon::Float64,
    number_options::Vector{Int64}
)::Matrix{Float64}
    s_final = zeros(length(P0_options), length(number_options))
    number_decision::Int64 = floor(time_horizon / decision_step)

    for (idx_P0, P0) in enumerate(P0_options)
        for (idx_number_option, number_option) in enumerate(number_options)
            s_final[idx_P0, idx_number_option] = _entropy(P0, influx, decision_step, number_decision, true, number_option)
        end
    end
    return s_final
end
