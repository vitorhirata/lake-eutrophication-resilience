function run_entropy(
    P0::Float64,
    influx::Float64,
    decision_step_options::Vector{Float64},
    number_decision_options::Vector{Int64},
)::Matrix{Float64}
    s_final = zeros(length(decision_step_options), length(number_decision_options))

    for (idx_decision_step, decision_step) in enumerate(decision_step_options)
        for (idx_number_decision, number_decision) in enumerate(number_decision_options)
            s_final[idx_decision_step, idx_number_decision] = _entropy(P0, influx, decision_step, number_decision)
        end
    end
    return s_final
end

function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    number_decision_options::Vector{Int64},
)::Matrix{Float64}
    s_final = zeros(length(P0_options), length(number_decision_options))

    for (idx_P0, P0) in enumerate(P0_options)
        for (idx_number_decision, number_decision) in enumerate(number_decision_options)
            s_final[idx_P0, idx_number_decision] = _entropy(P0, influx, decision_step, number_decision)
        end
    end
    return s_final
end

function run_entropy(
    P0_options::Vector{Float64},
    influx::Float64,
    decision_step::Float64,
    number_decision::Int64,
    number_options::Vector{Int64}
)::Matrix{Float64}
    s_final = zeros(length(P0_options), length(number_options))

    for (idx_P0, P0) in enumerate(P0_options)
        for (idx_number_option, number_option) in enumerate(number_options)
            s_final[idx_P0, idx_number_option] = _entropy(P0, influx, decision_step, number_decision, number_option)
        end
    end
    return s_final
end
