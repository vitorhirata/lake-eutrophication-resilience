function entropy(
        P0::Float64,
        I::Float64,
        decision_step::Float64,
        number_decision::Int64,
        prob::Float64 = 1.0;
        deterministic::Bool = true,
        max_options::Int64 = 10,
        minimum_influx::Float64 = 0.04,
        maximum_influx::Float64 = 0.30,
        past_P::Float64 = P0,
        method::String = "equal_probability"
)::Float64
    possible_influx = _possible_influx(P0, max_options, minimum_influx, maximum_influx)
    step_prob = _influx_probability(possible_influx, I, method, P0-past_P)

    final_prob = prob * step_prob
    if number_decision == 1
        return mapreduce(_causal_entropy, +, final_prob)
    end

    P_final = map(new_I -> _evolve_step(P0, new_I, decision_step, deterministic), possible_influx)
    results = map(input -> entropy(input[1], input[2], decision_step, number_decision-1, input[3];
                                    deterministic=deterministic, max_options=max_options,
                                    minimum_influx=minimum_influx, maximum_influx=maximum_influx,
                                    past_P=P0, method=method),
                  zip(P_final, possible_influx, final_prob))

    return sum(results)
end

function number_possible_influx(
    P::Float64,
    max_options::Int64,
    P_threshold::Float64 = 0.4,
    options_in_threshold::Int64 = max_options,
    final_P::Float64 = 3.0
)::Int64

    if P < P_threshold
        result = max_options + ((options_in_threshold - max_options) / P_threshold) * P
    elseif P < final_P
        result = options_in_threshold + ((1 - options_in_threshold) / (final_P - P_threshold)) * (P - P_threshold)
    else
        result = 1
    end

    return round(Int64, result)
end

function _possible_influx(
        P::Float64, max_number_options::Int64, minimum_influx::Float64 = 0.04, maximum_influx::Float64 = 0.30,
)::Vector{Float64}
    total_possible_influx = range(minimum_influx, maximum_influx, max_number_options)

    return collect(total_possible_influx[1:number_possible_influx(P, max_number_options)])
end

function _influx_probability(
        possible_influx::Vector{Float64}, influx::Float64, method::String, state_change::Float64
)::Vector{Float64}
    if method == "equal_probability"
        step_prob = _influx_probability_simple(possible_influx)
    elseif method == "closer_more_likely"
        step_prob = _influx_probability_closer(possible_influx, influx)
    elseif method == "further_more_likely"
        step_prob = _influx_probability_further(possible_influx, influx)
    elseif method == "state_change"
        step_prob = _influx_probability_change(possible_influx, state_change)
    elseif method == "favor_positive"
        step_prob = _influx_probability_positive(possible_influx)
    elseif method == "favor_negative"
        step_prob = _influx_probability_negative(possible_influx)
    else
        error("invalid method in entropy")
    end
    return step_prob
end

function _causal_entropy(prob::Float64)::Float64
    if prob == 0
        return 0
    end
    return - (prob * log(prob))
end
