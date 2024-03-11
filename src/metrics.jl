function max_kendall_tau_idx(time_series::NamedDimsArray, times::StepRangeLen{Float64})::Tuple{Int64, Float64}
    kendall_range = _kendall_range(time_series, times)
    times = collect(times)
    time_series = parent(time_series)

    kendall = [corkendall(times[idx:end], time_series[idx:end]) for idx in kendall_range]
    max_kendall = findmax(abs.(kendall))
    return kendall_range[max_kendall[2]], max_kendall[1]
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

function _kendall_range(time_series::NamedDimsArray, times::StepRangeLen{Float64})::StepRangeLen{Int64}
    begin_kendall::Int64 = findfirst(x -> x != 0, time_series)
    step_kendall::Int64 = (2 / step(times))
    end_kendall::Int64 = length(time_series) - 2 * step_kendall
    return begin_kendall:step_kendall:end_kendall
end
