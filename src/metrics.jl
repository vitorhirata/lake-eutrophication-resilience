function finite_difference(time_series::NamedDimsArray, interval::Float64)::NamedDimsArray
    result = NamedDimsArray{(:P0, :time_horizon)}(zeros(size(time_series, :P0)-1, size(time_series, :time_horizon)))
    for horizon_idx in 1:size(time_series, :time_horizon)
        result[time_horizon=horizon_idx] = abs.(diff(parent(time_series[time_horizon=horizon_idx]))) / interval
    end
    return result
end

function threshold_points(p::NamedDimsArray, s::NamedDimsArray, times::StepRangeLen{Float64},
        variance::NamedDimsArray, autocorrelation::NamedDimsArray,
        base_influx::Float64, influx_tax::Float64, I_step::Float64 = 1.0
)::Tuple{Dict, Dict}
    points = Dict()
    kendall_tau = Dict()

    points[:p] = _cross_threshold(p, times, base_influx, influx_tax, I_step)
    points[:var], kendall_tau[:var] = _max_kendall_tau(variance, times)
    points[:autocorr], kendall_tau[:autocorr] = _max_kendall_tau(autocorrelation, times)

    points[:s], kendall_tau[:s] = zeros(Int64, size(s, :time_horizon)), zeros(size(s, :time_horizon))
    for idx in 1:size(s, :time_horizon)
        points[:s][idx], kendall_tau[:s][idx] = _max_kendall_tau(s[time_horizon=idx], times)
    end
    return points, kendall_tau
end

function normalize_pd(time_series::NamedDimsArray)::NamedDimsArray
    for index_time_horizon in 1:1:size(time_series, :time_horizon)
        time_series[time_horizon=index_time_horizon] /= maximum(time_series[time_horizon=index_time_horizon])
    end

    return time_series
end

function compute_variance(time_series::NamedDimsArray, variance_step::Int64)::NamedDimsArray
    variance_ts = NamedDimsArray{(:time, )}(zeros(length(time_series)))

    for index_time in (variance_step+1):1:length(time_series)
        variance_ts[index_time] = var(time_series[(index_time-variance_step):index_time])
    end

    return variance_ts
end

function compute_autocorrelation(time_series::NamedDimsArray, autocor_step::Int64)::NamedDimsArray
    autocor_ts = NamedDimsArray{(:time, )}(zeros(length(time_series)))

    for index_time in (autocor_step+1):1:length(time_series)
        autocor_ts[index_time] = cor(
            time_series[(index_time-autocor_step):(index_time-1)],
            time_series[(index_time-autocor_step+1):(index_time)])
    end

    return autocor_ts
end

function detrend(time_series::NamedDimsArray, times::StepRangeLen{Float64}, type::String)::NamedDimsArray
    if type == "linear"
        return linear_detrend(time_series, times)
    elseif type == "loess"
        return loess_detrend(time_series, times)
    else
        throw(ArgumentError("Invalid type of detrend"))
    end
end

function linear_detrend(time_series::NamedDimsArray, times::StepRangeLen{Float64})::NamedDimsArray
    result = NamedDimsArray{(:time, )}(zeros(length(time_series)))
    times = collect(times)

    for index_influx in 1:1:length(time_series)
        p_influx = parent(time_series[influx_tax=index_influx])
        data = (;times,p_influx)
        model = lm(@formula(p_influx ~ times), data)
        result[influx_tax=index_influx] = residuals(model)
    end
    return result
end

function loess_detrend(time_series::NamedDimsArray, times::StepRangeLen{Float64})::NamedDimsArray
    result = NamedDimsArray{(:time, )}(zeros(length(time_series)))
    times = collect(times)

    model = loess(times, parent(time_series), span=0.5)
    smooth_function = predict(model, times)
    result = time_series - smooth_function
    return result
end

function _cross_threshold(time_series::NamedDimsArray, times::StepRangeLen{Float64},
    base_influx::Float64, influx_tax::Float64, I_step::Float64 = 1.0
)
    influx_ts = [base_influx + influx_tax * (time ÷ I_step) for time in times]
    thresholds = map(influx -> get_root(1.3, influx), influx_ts)
    first_cross = findfirst(item -> item > 0, time_series - thresholds)
    return first_cross
end

function _max_kendall_tau(
    time_series::NamedDimsArray, times::StepRangeLen{Float64}, kendall_threshold::Float64 = 0.55
)::Tuple{Int64, Float64}
    begin_point, kendall_range = _kendall_range(time_series, times)
    times = collect(times)
    time_series = parent(time_series)

    kendall = [corkendall(times[begin_point:idx], time_series[begin_point:idx]) for idx in kendall_range]
    idx_threshold = findfirst(x -> x > kendall_threshold, abs.(kendall))

    if idx_threshold == nothing
        return length(times), round(kendall[end]; digits=2)
    end
    return kendall_range[idx_threshold], round(kendall[idx_threshold]; digits=2)
end
#    max_kendall = findmax(abs.(kendall))
#    return kendall_range[max_kendall[2]], round(max_kendall[1]; digits=2)

function _kendall_range(
    time_series::NamedDimsArray, times::StepRangeLen{Float64}, kendell_shift::Int64 = 20
)::Tuple{Int64, StepRangeLen{Int64}}
    step_kendall::Int64 = (1 / step(times))
    begin_point::Int64 = findfirst(x -> x != 0, time_series)
    begin_kendall::Int64 = begin_point + kendell_shift * step_kendall
    end_kendall::Int64 = length(time_series) - 2 * step_kendall
    return begin_point, begin_kendall:step_kendall:end_kendall
end
