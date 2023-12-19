function compute_variance(time_series::Vector{Any}, variance_step::Int64)::Matrix{Float64}
    variance_ts = zeros(length(time_series[1]), length(time_series))

    for index_influx in 1:1:length(time_series)
        for index_time in (variance_step+1):1:length(time_series[1])
            variance_ts[index_time, index_influx] = var(time_series[index_influx][(index_time-variance_step):index_time])
        end
    end

    return variance_ts
end

function compute_autocorrelation(time_series::Vector{Any}, autocor_step::Int64)::Matrix{Float64}
    autocor_ts = zeros(length(time_series[1]), length(time_series))

    for index_influx in 1:1:length(time_series)
        for index_time in (autocor_step+1):1:length(time_series[1])
            autocor_ts[index_time, index_influx] = autocor(time_series[index_influx][(index_time-autocor_step):index_time], [1])[1]
        end
    end

    return autocor_ts
end
