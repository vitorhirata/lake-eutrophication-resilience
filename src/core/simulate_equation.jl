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

function get_root(root, parameter)
    Z = ZeroProblem(_f_root, root)
    sol = solve(Z, Order1(), p=parameter)
    return sol
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
