function compute_number_decision(time_horizon::Float64, decision_steps::Vector{Float64})
    result::Vector{Int64} = compute_number_decision.(time_horizon, decision_steps)
    return result
end

function compute_number_decision(time_horizons::Vector{Float64}, decision_step::Float64)
    result::Vector{Int64} = compute_number_decision.(time_horizons, decision_step)
    return result
end

function compute_number_decision(time_horizon::Float64, decision_step::Float64)::Int64
    result::Int64 = floor(time_horizon / decision_step)
    return result
end
