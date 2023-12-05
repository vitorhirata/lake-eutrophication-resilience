using DifferentialEquations

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

function number_possible_influx(
    P::Float64,
    P_threshold::Float64 = 0.75,
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

function _entropy(P0::Float64, a::Float64, step::Float64, time_limit::Int64, prob::Float64 = 1.0)::Float64
    possible_a_vec = _possible_influx(P0)
    step_prob = 1.0 / length(possible_a_vec)

    if time_limit == 1
        final_prob = prob * step_prob
        return length(possible_a_vec) * ( - final_prob * log(final_prob))
    end

    final_prob = prob * step_prob
    P_final = map(a -> _evolve_step(P0, a, step), possible_a_vec)
    results = map(input -> _entropy(input[1], input[2], step, time_limit-1, final_prob), zip(P_final, possible_a_vec))

    return sum(results)
end

function _evolve_step(P0::Float64, a::Float64, step::Float64)::Float64
    tspan = (0.0, step)
    prob = ODEProblem(_f, P0, tspan, a)
    sol = solve(prob, Tsit5())
    return sol.u[end]
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

function _possible_influx(
    P::Float64,
    total_possible_influx::StepRangeLen = 0.0:0.04:0.36,
)::Vector{Float64}
    return collect(total_possible_influx[1:number_possible_influx(P)])
end

# Function not being used
function run_scenarios(
        P0::Float64, a::Float64, t_max::Float64, step::Float64, time_limit::Int64, n_scenarious::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}

    times = StepRange(1, step, t_max)
    s_final = zeros(length(times))
    P_final = zeros(length(times))

    for _ in 1:n_scenarious
        s, P = _run_scenario(P0, a, t_max, step, time_limit)
        s_final += s
        P_final += P
    end

    s_final /= n_scenarious
    P_final /= n_scenarious

    return s_final, y_final
end

# Function not being used
function _run_scenario(
        P0::Float64, a::Float64, t_max::Float64, step::Float64, time_limit::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}

    times = StepRange(1, step, t_max)
    s_final = zeros(length(times))
    P_final = zeros(length(times))

    for (index, t0) in enumerate(times)
        P_final[index] = P0
        s_final[index] = _entropy(P0, a, step, time_limit)

        if t0 != 1
            a = rand(_possible_influx(P0))
        end

        P0 = _evolve_step(P0, a, step)
    end

    return s_final, y_final
end
