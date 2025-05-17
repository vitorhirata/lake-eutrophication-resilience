function distance_basin_threshold(
        P0_options::Vector{Float64},
        influx::Float64,
        time_horizons::Vector{Float64},
        decision_step::Float64,
        timestamp::String,
        one_plot::Bool=true
)
    s = readdlm("../output/$(timestamp)_distance_basin_s.csv", ',')
    s = NamedDimsArray{(:P0, :time_horizon)}(s)
    number_decision = compute_number_decision(time_horizons, decision_step)
    s = PathwayDiversity.normalize_pd(s, number_decision)
    s_detrended = PathwayDiversity.detrend(s, P0_options)
    s_diff = PathwayDiversity.finite_difference(s_detrended, P0_options[2]-P0_options[1])
    peaks_idx = PathwayDiversity.find_peaks(s_diff, P0_options[2:end])

    _plot_distance_threshold(s, s_diff, P0_options, time_horizons, timestamp, peaks_idx, one_plot)
end

function states_distribution(P0_options::Vector{Float64}, time_horizon::Float64, decision_step::Float64,
        timestamp::String, influx::Float64 = 0.1
)
    number_decision = compute_number_decision(time_horizon, decision_step)
    _plot_states_distribution(P0_options, number_decision, timestamp)
end

function early_warning_signals(
        timestamp::String;
        P0::Float64,
        influx::Float64,
        influx_tax::Float64,
        times::StepRangeLen{Float64},
        decision_step::Float64,
        time_horizons::Vector{Float64}
)

    p_ts = readdlm("../output/$(timestamp)_ews_p.csv", ',')
    p_ts = NamedDimsArray{(:time,)}(vec(p_ts))
    s_ts = readdlm("../output/$(timestamp)_ews_s.csv", ',')
    s_ts = NamedDimsArray{(:time, :time_horizon)}(s_ts)

    residuals = PathwayDiversity.detrend(p_ts,times, "loess")
    number_decision = compute_number_decision(time_horizons, decision_step)
    s_ts = PathwayDiversity.normalize_pd(s_ts, number_decision)

    # Compute variance
    variance_time_step = 30
    variance_idx_step::Int64 = variance_time_step ÷ step(times)
    variance_ts = PathwayDiversity.compute_variance(residuals, variance_idx_step)

    # Compute autocorrelation
    autocorr_time_step = 30
    autocorr_idx_step::Int64 = autocorr_time_step ÷ step(times)
    autocorr_ts = PathwayDiversity.compute_autocorrelation(residuals, autocorr_idx_step)

    # Compute threshold
    threshold_idx = PathwayDiversity.cross_threshold(p_ts, times, influx, influx_tax)

    # Compute kendall_tau
    kendall_tau = PathwayDiversity.kendall_tau(s_ts, variance_ts, autocorr_ts, times, threshold_idx)

    _plot_early_warning_signals(timestamp, p_ts, s_ts, variance_ts, autocorr_ts, time_horizons, times,
                                threshold_idx, kendall_tau)
    _plot_early_warning_residuals(timestamp, p_ts, residuals, times, threshold_idx)
end

function sensitivity(
        P0_options::Vector{Float64},
        influx::Float64,
        time_horizon::Float64,
        scenarios::Vector{Dict{Symbol, Any}},
        timestamp::String,
        relative::Bool = false
    )
    s = readdlm("../output/$(timestamp)_sensitivity_s.csv", ',')
    s = NamedDimsArray{(:P0, :type)}(s)

    decision_step = broadcast(scenario -> scenario[:decision_step], scenarios)
    number_decision = compute_number_decision(time_horizon, decision_step)
    s = PathwayDiversity.normalize_pd(s, number_decision)
    if relative
        s = PathwayDiversity.relative_pd(s)
    end

    _plot_sensitivity(s, P0_options, timestamp, scenarios, relative)
end
