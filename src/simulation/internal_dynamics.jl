# Compute the evolution of the state probability distribution over time for different initial conditions
function run_state_distribution(
    P0_options::Vector{Float64},
    time_horizon::Float64,
    decision_step::Float64
)::String
    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_state_distribution_"

    for (idx_P0, P0) in enumerate(P0_options)
        states_distribution = _states_distribution(P0, time_horizon, decision_step)
        println("Finished model for P0=$(P0)")
        writedlm("$(base_filename)$(idx_P0).csv",  states_distribution, ',')
    end
    return timestamp
end

# Evolve one time step and compute the number of options available for each P0 and influx.
function run_number_options_simulation(
    P0_options::Vector{Float64},
    influx_options::Vector{Float64},
    max_number_options::Int64 = 20,
    decision_step::Float64 = 5.0,
    deterministic::Bool = true
)::NamedDimsArray
    result = NamedDimsArray{(:P0, :influx)}(zeros(length(P0_options), length(influx_options)))

    for (idx_P0, P0) in enumerate(P0_options), (idx_influx, influx) in enumerate(influx_options)
        new_P = _evolve_step(P0, influx, decision_step, deterministic)
        result[P0=idx_P0, influx=idx_influx] = number_possible_influx(new_P, max_number_options)
    end
    return result
end

# Compute the range of states reachable from initial states
function run_range_states_simulation(
    P0_options::Vector{Float64},
    max_number_options::Int64 = 10,
    decision_step::Float64 = 5.0,
    deterministic::Bool = true
)::NamedDimsArray
    result = NamedDimsArray{(:P0, :new_P)}(zeros(length(P0_options), 2))

    for (idx, P0) in enumerate(P0_options)
        possible_a_vec = _possible_influx(P0, max_number_options)

        lower_P = _evolve_step(P0, possible_a_vec[begin], decision_step, deterministic)
        higher_P = _evolve_step(P0, possible_a_vec[end], decision_step, deterministic)
        result[P0=idx, new_P=1] = lower_P
        result[P0=idx, new_P=2] = higher_P
    end
    return result
end

function _states_distribution(
    P0::Float64,
    time_horizon::Float64,
    decision_step::Float64,
    max_number_options::Int64 = 10,
    deterministic::Bool = true
)::Vector{Vector{Float64}}

    states_distribution = Vector{Vector{Float64}}()
    number_decision = compute_number_decision(time_horizon, decision_step)
    push!(states_distribution, [P0])

    for decision in 2:(number_decision+1)
        push!(states_distribution, [])
        for state in states_distribution[decision-1]
            possible_a_vec = _possible_influx(state, max_number_options)
            P_final = map(new_I -> _evolve_step(state, new_I, decision_step, deterministic), possible_a_vec)
            append!(states_distribution[decision], P_final)
        end
    end
    return states_distribution
end
