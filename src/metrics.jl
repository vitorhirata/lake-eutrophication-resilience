function compute_variance(time_series::Matrix{Float64}, variance_step::Int64)::Matrix{Float64}
    variance_ts = zeros(size(time_series))

    for index_influx in 1:1:size(time_series, 1), index_time in (variance_step+1):1:size(time_series, 2)
        variance_ts[index_influx, index_time] = var(time_series[index_influx, (index_time-variance_step):index_time])
    end

    return variance_ts
end

function compute_autocorrelation(time_series::Matrix{Float64}, autocor_step::Int64)::Matrix{Float64}
    autocor_ts = zeros(size(time_series))

    for index_influx in 1:1:size(time_series, 1), index_time in (autocor_step+1):1:size(time_series, 2)
        autocor_ts[index_influx, index_time] = cor(
            time_series[index_influx, (index_time-autocor_step):(index_time-1)],
            time_series[index_influx, (index_time-autocor_step+1):(index_time)])
    end

    return autocor_ts
end
