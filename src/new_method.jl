function _new_entropy(
        P0::Float64,
        influx::Float64,
        decision_step::Float64,
        number_decision::Int64;
        deterministic::Bool = true,
        max_options::Int64 = 10,
        minimum_influx::Float64 = 0.04,
        maximum_influx::Float64 = 0.30,
        method::String = "equal_probability",
        possible_states::StepRangeLen{Float64} = 0:0.01:4,
        matrices::Union{Nothing, Vector{Matrix{Float64}}} = nothing,
)::Float64

    possible_influx = range(minimum_influx, maximum_influx, max_options)
    if isnothing(matrices)
        matrices = _full_state_transition_matrix(possible_influx, possible_states, decision_step, deterministic, method)
    end
    past_state = _initial_state_matrix(P0, influx, possible_states, possible_influx)
    entropy = 0

    for decision in 1:number_decision
        transition_matrix = _option_transition_matrix(past_state, matrices)
        entropy += _increment_entropy(past_state, transition_matrix)
        past_state = PathwayDiversity._iterate_state_matrix(past_state, matrices)
    end

    return entropy
end

function _initial_state_matrix(
        P0::Float64, influx::Float64, possible_states::StepRangeLen{Float64}, possible_influx::StepRangeLen{Float64}
)::Matrix{Float64}
    matrix = zeros(length(possible_states), length(possible_influx))

    idx_state = searchsortedfirst(possible_states, P0, lt = <=) - 1
    idx_influx = searchsortedfirst(possible_influx, influx, lt = <=) - 1
    matrix[idx_state,idx_influx] = 1.0
    return matrix
end

function _option_transition_matrix(past_state::Matrix{Float64}, matrices::Vector{Matrix{Float64}})::Matrix{Float64}
    option_matrix = zeros(size(past_state, 2), size(past_state, 2))

    for row in 1:size(past_state, 2)
        option_matrix[row,:] = sum(matrices[row] * past_state, dims = 1)
    end

    normalisation = map(el -> iszero(el) ? 1 : 1 / el, sum(option_matrix, dims=1))
    for col in 1:length(normalisation)
        option_matrix[:,col] *= normalisation[col]
    end

    return option_matrix
end


function _iterate_state_matrix(past_state::Matrix{Float64}, matrices::Vector{Matrix{Float64}})::Matrix{Float64}
    new_state = zeros(size(past_state))
    current_x = sum(past_state, dims = 2)
    for col in 1:size(past_state, 2)
        new_state[:,col] = matrices[col] * current_x
    end

    return new_state
end

function _full_state_transition_matrix(
    possible_influx::StepRangeLen{Float64} = range(0.04, 0.30, 10), possible_states::StepRangeLen{Float64} = 0:0.05:4,
    decision_step::Float64 = 5.0, deterministic::Bool=true, method::String = "equal_probability"
)::Vector{Matrix{Float64}}
    matrices = map(possible_influx) do influx
        _state_transition_matrix(possible_states, influx, length(possible_influx), decision_step, deterministic, method)
    end
    return matrices
end

function _state_transition_matrix(
    possible_x::StepRangeLen{Float64}, influx::Float64, max_option::Int64, decision_step::Float64, deterministic::Bool,
    method::String
)::Matrix{Float64}
    matrix = zeros(length(possible_x), length(possible_x))
    for (idx_begin_state, begin_state) in enumerate(possible_x)
        possible_influx = PathwayDiversity._possible_influx(begin_state, max_option)

        !(influx in possible_influx) && continue

        probabilities = PathwayDiversity._influx_probability(possible_influx, influx, method, -1.0)

        idx_prob = findfirst(isequal(influx), possible_influx)
        probability = probabilities[idx_prob]
        end_state = PathwayDiversity._evolve_step(begin_state, influx, decision_step, deterministic)
        end_bin = searchsortedfirst(possible_x, end_state) - 1

        matrix[end_bin,idx_begin_state] += probability
    end
    return matrix
end

function _increment_entropy(past_state::Matrix{Float64}, transition_matrix::Matrix{Float64})::Float64
    option_diversity = sum(PathwayDiversity._causal_entropy.(transition_matrix), dims = 1)
    influx_state = sum(past_state, dims = 1)
    return dot(influx_state, option_diversity)
end
