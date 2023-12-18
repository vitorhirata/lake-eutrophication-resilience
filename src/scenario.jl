function run_scenarios(
        P0::Float64,
        I::Float64,
        t_max::Int64,
        step::Float64,
        decision_step::Float64,
        time_horizon::Float64,
        n_scenarious::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}

    times = 1:step:t_max
    s_final = zeros(length(times))
    P_final = zeros(length(times))
    number_decision::Int64 = floor(time_horizon / decision_step)

    for _ in 1:n_scenarious
        s, P = _run_scenario(P0, I, t_max, step, decision_step, number_decision)
        s_final += s
        P_final += P
    end

    s_final /= n_scenarious
    P_final /= n_scenarious

    return s_final, P_final
end

function _run_scenario(
        P0::Float64, I::Float64, t_max::Int64, step::Float64, decision_step::Float64, number_decision::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}

    times = 1:step:t_max
    s_final = zeros(length(times))
    P_final = zeros(length(times))

    for (index, t0) in enumerate(times)
        P_final[index] = P0
        s_final[index] = _entropy(P0, I, decision_step, number_decision, false)

        P0 = _evolve_step(P0, I, step, false)
    end

    return s_final, P_final
end
