function run_scenario(
        P0::Float64,
        I::Float64,
        I_tax::Float64,
        times::StepRangeLen{Float64},
        decision_step::Float64,
        time_horizons::Vector{Float64},
        I_step::Float64 = 1.0
)::Tuple{NamedDimsArray, NamedDimsArray}

    p = NamedDimsArray{(:time,)}(zeros(length(times)))
    s = NamedDimsArray{(:time, :time_horizon)}(zeros(length(times), length(time_horizons)))
    number_decision::Vector{Int64} = map(time_horizon -> floor(time_horizon / decision_step), time_horizons)

    for (index_t, t0) in enumerate(times)
        p[index_t] = P0

        for (index_h, time_horizon) in enumerate(time_horizons)
            s[index_t, index_h] = _entropy(P0, I, decision_step, number_decision[index_h], false)
        end

        P0 = _evolve_step(P0, I, step(times), false)
        if t0 % I_step == 0
            I += I_tax
        end
    end

    return p, s
end
