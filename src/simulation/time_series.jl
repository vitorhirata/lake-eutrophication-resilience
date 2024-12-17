function run_time_series(;
        P_init::Float64,
        influx::Float64,
        influx_tax::Float64,
        times::StepRangeLen{Float64},
        decision_step::Float64,
        time_horizons::Vector{Float64},
        influx_step::Float64 = 1.0
)::String
    p = NamedDimsArray{(:time,)}(zeros(length(times)))
    influx_ts = NamedDimsArray{(:time,)}(zeros(length(times)))
    s = NamedDimsArray{(:time, :time_horizon)}(zeros(length(times), length(time_horizons)))
    number_decision = compute_number_decision(time_horizons, decision_step)

    for (index_t, t0) in enumerate(times)
        p[index_t] = P_init
        influx_ts[index_t] = influx

        for (index_h, time_horizon) in enumerate(time_horizons)
            s[index_t, index_h] = entropy(P_init, influx, decision_step, number_decision[index_h]; deterministic=false)
        end

        P_init = _evolve_step(P_init, influx, step(times), false)
        if t0 % influx_step == 0
            influx += influx_tax
        end
        t0 % 20 == 0 && println("Finished running $(t0) time step")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_ews_"
    writedlm("$(base_filename)p.csv",  p, ',')
    writedlm("$(base_filename)influx.csv",  influx_ts, ',')
    writedlm("$(base_filename)s.csv",  s, ',')
    return timestamp
end
