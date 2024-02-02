function run_scenarios(
        P0::Float64,
        I::Float64,
        I_taxes::Vector{Float64},
        t_max::Int64,
        step::Float64,
        decision_step::Float64,
        time_horizon::Float64,
)::NamedDimsArray

    times = 1:step:t_max
    results = NamedDimsArray{(:time, :influx_tax, :type)}(zeros(length(times), length(I_taxes), 2))
    number_decision::Int64 = floor(time_horizon / decision_step)

    for (index, I_tax) in enumerate(I_taxes)
        results[influx_tax=index] = _run_scenario(P0, I, I_tax, number_decision, t_max, step, decision_step)
    end

    return results
end

function _run_scenario(
        P0::Float64, I::Float64, I_tax::Float64, number_decision::Int64,
        t_max::Int64, step::Float64, decision_step::Float64, I_step::Float64 = 1.0
)::NamedDimsArray

    times = 1:step:t_max
    result = NamedDimsArray{(:time, :type)}(zeros(length(times), 2))

    for (index, t0) in enumerate(times)
        result[index, 1] = P0
        result[index, 2] = _entropy(P0, I, decision_step, number_decision, false)

        P0 = _evolve_step(P0, I, step, false)
        if t0 % I_step == 0
            I += I_tax
        end
    end

    return result
end
