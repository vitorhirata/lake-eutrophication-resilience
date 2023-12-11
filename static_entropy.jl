
function run_entropy(
    P0::Float64,
    influx::Float64,
    step_options::Vector{Float64},
    time_limit_options::Vector{Int64},
)::Matrix{Float64}
    s_final = zeros(length(step_options), length(time_limit_options))

    for (idx_step, step) in enumerate(step_options)
        for (idx_time, time_limit) in enumerate(time_limit_options)
            s_final[idx_step, idx_time] = _entropy(P0, influx, step, time_limit)
        end
    end
    return s_final
end

function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    step::Float64,
    time_limit_options::Vector{Int64},
)::Matrix{Float64}
    s_final = zeros(length(P0_options), length(time_limit_options))

    for (idx_P0, P0) in enumerate(P0_options)
        for (idx_time, time_limit) in enumerate(time_limit_options)
            s_final[idx_P0, idx_time] = _entropy(P0, influx, step, time_limit)
        end
    end
    return s_final
end

