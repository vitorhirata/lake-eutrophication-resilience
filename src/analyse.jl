function early_warning_signals(
        timestamp::String;
        P_init::Float64,
        influx::Float64,
        influx_tax::Float64,
        times::StepRangeLen{Float64},
        decision_step::Float64,
        time_horizons::Vector{Float64}
)

    p = readdlm("../output/$(timestamp)_ews_p.csv", ',')
    p = NamedDimsArray{(:time,)}(vec(p))
    s = readdlm("../output/$(timestamp)_ews_s.csv", ',')
    s = NamedDimsArray{(:time, :time_horizon)}(s)

    residuals = PathwayDiversity.detrend(p, times, "loess")
    s = PathwayDiversity.normalize_pd(s)

    # Compute variance
    variance_time_step = 20
    variance_idx_step::Int64 = variance_time_step ÷ step(times)
    variance_ts = PathwayDiversity.compute_variance(residuals, variance_idx_step)

    # Compute autocorrelation
    autocorr_time_step = 20
    autocorr_idx_step::Int64 = autocorr_time_step ÷ step(times)
    autocorr_ts = PathwayDiversity.compute_autocorrelation(residuals, autocorr_idx_step)

    tipping_points, kendall_tau = PathwayDiversity.threshold_points(p, s, times, variance_ts, autocorr_ts,
                                                                    influx, influx_tax)

    _plot_early_warning_signals(p, s, residuals, variance_ts, autocorr_ts, time_horizons,
                                times, variance_time_step, autocorr_time_step, tipping_points, kendall_tau)
end
