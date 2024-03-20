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

function early_warning_kendall_taus(
        timestamps::Vector{String};
        P_init::Float64,
        influx::Float64,
        influx_tax::Float64,
        times::StepRangeLen{Float64},
        decision_step::Float64,
        time_horizons::Vector{Float64}
)
    kendall_taus = NamedDimsArray{(:simulation, :data, :type)}(zeros(length(timestamps), 2+length(time_horizons), 2))

    for (idx, timestamp) in enumerate(timestamps)
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

        if tipping_points[:p] == nothing
            tipping_points[:p] = length(times)
        end
        kendall_taus[simulation=idx, data=1, type=1] = tipping_points[:var] / tipping_points[:p]
        kendall_taus[simulation=idx, data=1, type=2] = kendall_tau[:var]
        kendall_taus[simulation=idx, data=2, type=1] = tipping_points[:autocorr] / tipping_points[:p]
        kendall_taus[simulation=idx, data=2, type=2] = kendall_tau[:autocorr]
        for horizon in 1:length(time_horizons)
            kendall_taus[simulation=idx, data=2+horizon, type=1] = tipping_points[:s][horizon] / tipping_points[:p]
            kendall_taus[simulation=idx, data=2+horizon, type=2] = kendall_tau[:s][horizon]
        end
    end
    kendall_taus[type=2] *= step(times)

    _plot_kendall_taus(kendall_taus, time_horizons)
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

    _plot_distance_threshold(s, s_diff, distance_threshold, time_horizons)
end

function decision_scales(decision_steps::Vector{Float64}, time_horizons::Vector{Float64}, timestamp::String)
    s = readdlm("../output/$(timestamp)_decision_scale_s.csv", ',')
    s = NamedDimsArray{(:decision_step, :time_horizon)}(s)
    _plot_decision_scales(s, time_horizons, decision_steps)
end

function scaling(P_init_options::Vector{Float64}, number_options::Vector{Int64}, timestamp::String)
    s = readdlm("../output/$(timestamp)_scale_initial_s.csv", ',')
    s = NamedDimsArray{(:number_options, :P0)}(s)
    _plot_scaling(s, P_init_options, number_options)
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
