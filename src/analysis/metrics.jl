function finite_difference(time_series::NamedDimsArray, interval::Float64)::NamedDimsArray
    result = NamedDimsArray{(:P0, :time_horizon)}(zeros(size(time_series, :P0)-1, size(time_series, :time_horizon)))
    for horizon_idx in 1:size(time_series, :time_horizon)
        result[time_horizon=horizon_idx] = abs.(diff(parent(time_series[time_horizon=horizon_idx]))) / interval
    end
    return result
end

function kendall_tau(s::NamedDimsArray, variance::NamedDimsArray, autocorrelation::NamedDimsArray,
        times::StepRangeLen{Float64}, threshold_idx::Int64
)::Dict{Symbol, Float64}
    result = Dict{Symbol, Float64}()

    result[:var] = kendall_tau(variance, times, threshold_idx)
    result[:autocorr] = kendall_tau(autocorrelation, times, threshold_idx)
    for idx in 1:size(s, :time_horizon)
        result[Symbol(:s_, idx)] = kendall_tau(s[time_horizon=idx], times, threshold_idx)
    end
    return result
end

function kendall_tau(time_series::NamedDimsArray, times::StepRangeLen{Float64},
                     threshold_idx::Int64, kendall_tau_time_offset = 6
)::Float64
    delta_kendall_tau_idx = length(times) - threshold_idx + Int64(kendall_tau_time_offset / step(times))
    initial_time_idx = length(times) - length(time_series) + 1
    analysis_times = collect(times[initial_time_idx:(length(times) - delta_kendall_tau_idx)])
    analysis_time_series = parent(time_series[1:length(time_series) - delta_kendall_tau_idx])

    result = corkendall(analysis_times, analysis_time_series)
    return round(result; digits=2)
end

function find_peaks(time_series::NamedDimsArray, time::Vector{Float64})::Vector{Int64}
    peaks_idx = zeros(Int64, size(time_series, :time_horizon))
    if size(time_series, :time_horizon) == 1
        peak_prominance = [0.01] # [0.001, 0.01, 0.01, 0.01]
    else
        peak_prominance = [0.001, 0.01, 0.01, 0.01]
    end

    for horizon_idx in 1:size(time_series, :time_horizon)
        pks = findmaxima(time_series[time_horizon=horizon_idx]) |> peakproms!(; min=peak_prominance[horizon_idx])
        peaks_idx[horizon_idx] = pks[:indices][1]
    end
    return peaks_idx
end

function normalize_pd(time_series::NamedDimsArray, dim::Symbol)::NamedDimsArray
    for index_dim in 1:1:size(time_series, dim)
        eval(:( $time_series[$dim = $index_dim] /= maximum($time_series[$dim = $index_dim]) ))
    end

    return time_series
end

function normalize_pd(
        time_series::NamedDimsArray, number_decision::Vector{Int64}, max_options::Int64 = 10
)::NamedDimsArray
    for index in 1:1:size(time_series, 2)
        time_series[:,index] /= (number_decision[index] * max_options)
    end

    return time_series
end

function relative_pd(time_series::NamedDimsArray)::NamedDimsArray
    for index_dim in 2:1:size(time_series, :type)
        time_series[type = index_dim] = time_series[type = index_dim] ./ time_series[type = 1]
    end

    return time_series
end

function compute_variance(time_series::NamedDimsArray, variance_step::Int64)::NamedDimsArray
    variance_ts = NamedDimsArray{(:time, )}(zeros(length(time_series)))

    for index_time in (variance_step+1):1:length(time_series)
        variance_ts[index_time] = var(time_series[(index_time-variance_step):index_time])
    end

    return variance_ts[(variance_step+1):end]
end

function compute_autocorrelation(time_series::NamedDimsArray, autocor_step::Int64)::NamedDimsArray
    autocor_ts = NamedDimsArray{(:time, )}(zeros(length(time_series)))

    for index_time in (autocor_step+1):1:length(time_series)
        autocor_ts[index_time] = cor(
            time_series[(index_time-autocor_step):(index_time-1)],
            time_series[(index_time-autocor_step+1):(index_time)])
    end

    return autocor_ts[(autocor_step+1):end]
end

function detrend(time_series::NamedDimsArray, times::StepRangeLen{Float64}, type::String)::NamedDimsArray
    if type == "linear"
        return linear_detrend(time_series, times)
    elseif type == "loess"
        return loess_detrend(time_series, collect(times), :time, true)
    else
        throw(ArgumentError("Invalid type of detrend"))
    end
end

function detrend(time_series::NamedDimsArray, times::Vector{Float64})::NamedDimsArray
    result = NamedDimsArray{(:P0, :time_horizon)}(zeros(size(time_series, :P0), size(time_series, :time_horizon)))

    for idx_horizon in 1:size(time_series, :time_horizon)
        result[time_horizon=idx_horizon] = loess_detrend(time_series[time_horizon=idx_horizon], times, :P0, false)
    end
    return result
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

function loess_detrend(time_series::NamedDimsArray, times::Vector{Float64}, dim::Symbol, residual::Bool)::NamedDimsArray
    result = NamedDimsArray{(dim, )}(zeros(length(time_series)))

    model = loess(times, parent(time_series), span=0.5)
    smooth_function = predict(model, times)
    result = residual ? time_series - smooth_function : NamedDimsArray{(dim, )}(smooth_function)
    return result
end

function cross_threshold(time_series::NamedDimsArray, times::StepRangeLen{Float64},
    base_influx::Float64, influx_tax::Float64, I_step::Float64 = 1.0
)
    influx_ts = [base_influx + influx_tax * (time ÷ I_step) for time in times]
    thresholds = map(influx -> get_root(1.3, influx), influx_ts)
    first_cross = findfirst(item -> item > 0, time_series - thresholds)
    return first_cross
end

