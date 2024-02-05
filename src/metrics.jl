function compute_variance(time_series::NamedDimsArray, variance_step::Int64)::NamedDimsArray
    variance_ts = NamedDimsArray{(:time, :influx_tax)}(zeros(size(time_series)))

    for index_influx in 1:1:size(time_series, :influx_tax), index_time in (variance_step+1):1:size(time_series, :time)
        variance_ts[time=index_time, influx_tax=index_influx] =
            var(time_series[time=(index_time-variance_step):index_time, influx_tax=index_influx])
    end

    return variance_ts
end

function compute_autocorrelation(time_series::NamedDimsArray, autocor_step::Int64)::NamedDimsArray
    autocor_ts = NamedDimsArray{(:time, :influx_tax)}(zeros(size(time_series)))

    for index_influx in 1:1:size(time_series, :influx_tax), index_time in (autocor_step+1):1:size(time_series, :time)
        autocor_ts[time=index_time, influx_tax=index_influx] = cor(
            time_series[time=(index_time-autocor_step):(index_time-1), influx_tax=index_influx],
            time_series[time=(index_time-autocor_step+1):(index_time), influx_tax=index_influx])
    end

    return autocor_ts
end

function detrend(time_series::NamedDimsArray, times::Vector{Float64})
    result = NamedDimsArray{(:time, :influx_tax)}(zeros(size(time_series)))

    for index_influx in 1:1:size(time_series, :influx_tax)
        p_influx = parent(time_series[influx_tax=index_influx])
        data = (;times,p_influx)
        model = lm(@formula(p_influx ~ times), data)
        result[influx_tax=index_influx] = residuals(model)
    end
    return result
end
