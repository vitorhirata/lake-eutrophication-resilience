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

function get_root(root, parameter)
    Z = ZeroProblem(_f_root, root)
    sol = solve(Z, Order1(), p=parameter)
    return sol
end

function _entropy(
        P0::Float64,
        I::Float64,
        decision_step::Float64,
        number_decision::Int64,
        prob::Float64 = 1.0;
        deterministic::Bool = true,
        max_options::Int64 = 10
)::Float64
    possible_a_vec = _possible_influx(P0, max_options)
    step_prob = _influx_probability(possible_a_vec)

    final_prob = prob * step_prob
    if number_decision == 1
        return mapreduce(_causal_entropy, +, final_prob)
    end

    P_final = map(new_I -> _evolve_step(P0, new_I, decision_step, deterministic), possible_a_vec)
    results = map(input -> _entropy(input[1], input[2], decision_step, number_decision-1, input[3];
                                    deterministic=deterministic, max_options=max_options),
                  zip(P_final, possible_a_vec, final_prob))

    return sum(results)
end

function _possible_influx(P::Float64, max_number_options::Int64, maximun_influx::Float64 = 0.36)::Vector{Float64}
    total_possible_influx = range(0.0, maximun_influx, max_number_options)
    return collect(total_possible_influx[1:number_possible_influx(P, max_number_options)])
end

function _influx_probability(possible_influx::Vector{Float64})::Vector{Float64}
    prob = 1.0 / length(possible_influx)
    return fill(prob, length(possible_influx))
end


function _evolve_step(P0::Float64, I::Float64, step::Float64, deterministic::Bool)::Float64
    if deterministic
        _evolve_step_deterministic(P0, I, step)
    else
        _evolve_step_stochastic(P0, I, step)
    end
end

function _evolve_step_deterministic(P0::Float64, I::Float64, step::Float64)::Float64
    tspan = (0.0, step)
    prob = ODEProblem(_f, P0, tspan, I)
    sol = solve(prob, Tsit5())
    return sol.u[end]
end

function _evolve_step_stochastic(
        P0::Float64,
        I::Float64,
        step::Float64,
        μ::Float64 = 0.1,
        σ::Float64 = 0.01,
        dt::Float64 = 0.05,
)::Float64
    tspan = (0.0, step)
    W = GeometricBrownianMotionProcess(μ, σ, 0.0, 1.0, 1.0)
    prob = SDEProblem(_f, _g, P0, tspan, I, noise=W)
    sol = solve(prob, EM(), dt = dt)
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

_f_root(P, I) = _f(P, I, 0.0)
_g(P, I=nothing, time=nothing) = 1
_causal_entropy(prob) = - (prob * log(prob))
