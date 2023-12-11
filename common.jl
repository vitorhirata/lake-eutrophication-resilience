function number_possible_influx(
    P::Float64,
    P_threshold::Float64 = 0.6,
    possibilities_threshold::Int64 = 10,
    max_possibilities::Int64 = 10,
    final_P::Float64 = 3.0
)::Int64

    if P < P_threshold
        result = max_possibilities + ((possibilities_threshold - max_possibilities) / P_threshold) * P
    elseif P < final_P
        result = possibilities_threshold + ((1 - possibilities_threshold) / (final_P - P_threshold)) * (P - P_threshold)
    else
        result = 1
    end

    return round(Int64, result)
end

function _f(
    P::Float64,
    influx::Float64,
    time::Float64,
    # s::Float64 = 0.65,
    # e::Float64 = 2.5,
    # h::Float64 = 1.95,
)::Float64
    return influx - 0.65 * P + 2.5 * (P^2) / ((1.95)^2 + P^2)
end

function _evolve_step(P0::Float64, I::Float64, step::Float64)::Float64
    tspan = (0.0, step)
    prob = ODEProblem(_f, P0, tspan, I)
    sol = solve(prob, Tsit5())
    return sol.u[end]
end

function _entropy(P0::Float64, I::Float64, step::Float64, time_limit::Int64, prob::Float64 = 1.0)::Float64
    possible_a_vec = _possible_influx(P0)
    step_prob = 1.0 / length(possible_a_vec)

    if time_limit == 1
        final_prob = prob * step_prob
        return length(possible_a_vec) * ( - final_prob * log(final_prob))
    end

    final_prob = prob * step_prob
    P_final = map(a -> _evolve_step(P0, I, step), possible_a_vec)
    results = map(input -> _entropy(input[1], input[2], step, time_limit-1, final_prob), zip(P_final, possible_a_vec))

    return sum(results)
end

function _possible_influx(
    P::Float64,
    total_possible_influx::StepRangeLen = 0.0:0.04:0.36,
)::Vector{Float64}
    return collect(total_possible_influx[1:number_possible_influx(P)])
end
