function early_warning_signals(
        timestamp::String;
        P_init::Float64,
        influx::Float64,
        influx_tax::Float64,
        times::StepRangeLen{Float64},
        decision_step::Float64,
        time_horizons::Vector{Float64}
)

    p_ts = readdlm("../output/$(timestamp)_ews_p.csv", ',')
    p_ts = NamedDimsArray{(:time,)}(vec(p_ts))
    influx_ts = readdlm("../output/$(timestamp)_ews_influx.csv", ',')
    influx_ts = NamedDimsArray{(:time,)}(vec(influx_ts))
    s_ts = readdlm("../output/$(timestamp)_ews_s.csv", ',')
    s_ts = NamedDimsArray{(:time, :time_horizon)}(s_ts)

    residuals = PathwayDiversity.detrend(p_ts,times, "loess")
    s_ts = PathwayDiversity.normalize_pd(s_ts)

    # Compute variance
    variance_time_step = 30
    variance_idx_step::Int64 = variance_time_step ÷ step(times)
    variance_ts = PathwayDiversity.compute_variance(residuals, variance_idx_step)

    # Compute autocorrelation
    autocorr_time_step = 30
    autocorr_idx_step::Int64 = autocorr_time_step ÷ step(times)
    autocorr_ts = PathwayDiversity.compute_autocorrelation(residuals, autocorr_idx_step)

    # Compute distance to threshold
    thresholds = NamedDimsArray{(:time,)}(zeros(length(p_ts)))
    for (idx, infl) in enumerate(influx_ts)
        thresholds[idx] = p_ts[idx] - PathwayDiversity.get_root(1.3, infl)
    end

    tipping_points, kendall_tau = PathwayDiversity.threshold_points(p_ts,s_ts, times, variance_ts, autocorr_ts,
                                                                    influx, influx_tax)

    _plot_early_warning_signals(timestamp, p_ts,s_ts, residuals, variance_ts, autocorr_ts, influx_ts, thresholds,
                                time_horizons, times, variance_time_step, autocorr_time_step, tipping_points,
                                kendall_tau)
end

function distance_basin_threshold(
        P_init_options::Vector{Float64},
        influx::Float64,
        time_horizons::Vector{Float64},
        timestamp::String
)
    s = readdlm("../output/$(timestamp)_distance_basin_s.csv", ',')
    s = NamedDimsArray{(:P0, :time_horizon)}(s)
    s = PathwayDiversity.normalize_pd(s)
    s_diff = PathwayDiversity.finite_difference(s, 0.05)

    threshold = PathwayDiversity.get_root(1.3, influx)
    distance_threshold = threshold .- P_init_options
    peaks_idx = PathwayDiversity.find_peaks(s_diff, distance_threshold[2:end])

    _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons, timestamp, peaks_idx)
end

function decision_scales(decision_steps::Vector{Float64}, time_horizons::Vector{Float64}, timestamp::String)
    s = readdlm("../output/$(timestamp)_decision_scale_s.csv", ',')
    s = NamedDimsArray{(:decision_step, :time_horizon)}(s)
    _plot_decision_scales(s, time_horizons, decision_steps, timestamp)
end

function scaling(P_init_options::Vector{Float64}, number_options::Vector{Int64}, timestamp::String)
    s = readdlm("../output/$(timestamp)_scale_initial_s.csv", ',')
    s = NamedDimsArray{(:number_options, :P0)}(s)
    _plot_scaling(s, P_init_options, number_options, timestamp)
end

function number_options_simulation(
    P0_options::Vector{Float64}, influx_options::Vector{Float64}, max_number_options::Int64
)
    n_options = PathwayDiversity.run_number_options_simulation(P0_options, influx_options, max_number_options)
    _plot_number_options_simulation(P0_options, influx_options, n_options)
end

function range_states_simulation(P0_options::Vector{Float64}, max_number_options::Int64)
    range_states = PathwayDiversity.run_range_states_simulation(P0_options, max_number_options)
    _plot_range_states(P0_options, range_states)
end

function number_options(P_options::StepRangeLen{Float64}, max_number_options::Int64)
    n_possible_influx = PathwayDiversity.number_possible_influx.(P_options, max_number_options)
    _plot_number_option(n_possible_influx, P_options)
end

function bifurcation(influx_options_root::Vector{Float64})
    roots = zeros(3, length(influx_options_root))

    for (index, influx) in enumerate(influx_options_root)
        roots[1, index] = PathwayDiversity.get_root(0.1, influx)
        roots[2, index] = PathwayDiversity.get_root(1.3, influx)
        roots[3, index] = PathwayDiversity.get_root(2.7, influx)
    end
    _plot_bifurcation(roots, influx_options_root)
end
