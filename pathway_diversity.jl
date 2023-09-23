using DifferentialEquations

function run_entropy(
    x0::Float64,
    influx::Float64,
    step_options::Vector{Float64},
    time_limit_options::Vector{Int64},
)::Matrix{Float64}
    s_final = zeros(length(step_options), length(time_limit_options))

    for (idx_step, step) in enumerate(step_options)
        for (idx_time, time_limit) in enumerate(time_limit_options)
            s_final[idx_step, idx_time] = entropy(x0, influx, step, time_limit)
        end
    end
    return s_final
end

function run_entropy(
    x0_options::Vector{Float64},
    influx::Float64,
    step::Float64,
    time_limit_options::Vector{Int64},
)::Matrix{Float64}
    s_final = zeros(length(x0_options), length(time_limit_options))

    for (idx_x0, x0) in enumerate(x0_options)
        for (idx_time, time_limit) in enumerate(time_limit_options)
            s_final[idx_x0, idx_time] = entropy(x0, influx, step, time_limit)
        end
    end
    return s_final
end


function entropy(x0::Float64, a::Float64, step::Float64, time_limit::Int64, prob::Float64 = 1.0)::Float64
    possible_a_vec = _possible_a(x0)
    step_prob = 1.0 / length(possible_a_vec)

    if time_limit == 1
        final_prob = prob * step_prob
        return length(possible_a_vec) * ( - final_prob * log(final_prob))
    end

    final_prob = prob * step_prob
    x_final = map(a -> _evolve_step(x0, a, step), possible_a_vec)
    results = map(input -> entropy(input[1], input[2], step, time_limit-1, final_prob), zip(x_final, possible_a_vec))

    return sum(results)
end

function entropy_iterative()
    return 1
end

function run_scenarios(
        x0::Float64, a::Float64, t_max::Float64, step::Float64, time_limit::Int64, n_scenarious::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}

    times = StepRange(1, step, t_max)
    s_final = zeros(length(times))
    x_final = zeros(length(times))

    for _ in 1:n_scenarious
        s, x = _run_scenario(x0, a, t_max, step, time_limit)
        s_final += s
        x_final += x
    end

    s_final /= n_scenarious
    x_final /= n_scenarious

    return s_final, y_final
end

function _evolve_step(x0::Float64, a::Float64, step::Float64)::Float64
    tspan = (0.0, step)
    prob = ODEProblem(_f, x0, tspan, a)
    sol = solve(prob, Tsit5())
    return sol.u[end]
end

function _run_scenario(
        x0::Float64, a::Float64, t_max::Float64, step::Float64, time_limit::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}

    times = StepRange(1, step, t_max)
    s_final = zeros(length(times))
    x_final = zeros(length(times))

    for (index, t0) in enumerate(times)
        x_final[index] = x0
        s_final[index] = entropy(x0, a, step, time_limit)

        if t0 != 1
            a = rand(_possible_a(x0))
        end

        x0 = _evolve_step(x0, a, step)
    end

    return s_final, y_final
end

function _f(
    x::Float64,
    influx::Float64,
    time::Float64,
    # removal_rate::Float64 = 0.65,
    # recycling_rate::Float64 = 2.5,
    # recycling_half_saturation::Float64 = 1.95,
)::Float64
    return influx - 0.65 * x + 2.5 * (x^2) / ((1.95)^2 + x^2)
end

function _possible_a(
    x::Float64,
    x_limits::Tuple{Vararg{Float64}} = (1.5, 1.0, 0.5, -0.1),
    total_possible_a::Tuple{Vararg{Float64}} = (0.0, 0.1, 0.2, 0.3),
)::Tuple{Vararg{Float64}}

    number_a = _first_bigger_element_index(x, x_limits)
    return total_possible_a[1:number_a]
end

function _first_bigger_element_index(x::Float64, x_limits::Tuple{Vararg{Float64}})::Int64
    for (index, limit) in enumerate(x_limits)
        if x >= limit
            return index
        end
    end
end

