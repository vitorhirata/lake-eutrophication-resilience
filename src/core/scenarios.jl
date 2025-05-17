function _influx_probability_simple(possible_influx::Vector{Float64})::Vector{Float64}
    prob = 1.0 / length(possible_influx)
    return fill(prob, length(possible_influx))
end

function _influx_probability_closer(possible_influx::Vector{Float64}, past_influx::Float64)::Vector{Float64}
    past_influx_idx = findfirst(x -> past_influx <= x, possible_influx)
    if past_influx_idx == nothing
        past_influx_idx = length(possible_influx)
    end

    result = zeros(length(possible_influx))
    for idx in 1:length(possible_influx)
        result[idx] = _weight_probability_closer(idx, past_influx_idx)
    end

    return normalize(result, 1)
end

function _influx_probability_further(possible_influx::Vector{Float64}, past_influx::Float64)::Vector{Float64}
    past_influx_idx = findfirst(x -> past_influx <= x, possible_influx)
    if past_influx_idx == nothing
        past_influx_idx = length(possible_influx)
    end
    max_difference = maximum(abs.(collect(1:length(possible_influx)) .- past_influx_idx))

    result = zeros(length(possible_influx))
    for idx in 1:length(possible_influx)
        result[idx] = _weight_probability_further(idx, past_influx_idx, max_difference)
    end

    return normalize(result, 1)
end

function _influx_probability_positive(possible_influx::Vector{Float64})::Vector{Float64}
    result = reverse(range(1, step=4.0, length=length(possible_influx)))

    return normalize(collect(result), 1)
end

function _influx_probability_negative(possible_influx::Vector{Float64})::Vector{Float64}
    result = range(1, step=4.0, length=length(possible_influx))

    return normalize(collect(result), 1)
end

function _influx_probability_change(possible_influx::Vector{Float64}, state_change::Float64)::Vector{Float64}
    result::Vector{Float64} = []
    if state_change < -0.2
        result = range(1, step=4.0, length=length(possible_influx))
    elseif state_change < 0.2
        result = range(1, step=0.0, length=length(possible_influx))
    elseif state_change < 0.8
        result = reverse(range(1, step=4.0, length=length(possible_influx)))
    else
        result = reverse(range(1, step=10.0, length=length(possible_influx)))
    end
    return normalize(collect(result), 1)
end

function _weight_probability_closer(idx::Int64, reference::Int64)::Int64
    difference = abs(idx - reference)
    if difference == 0
        return 13
    elseif difference == 1
        return 9
    elseif difference == 2
        return 5
    elseif difference == 3
        return 1
    else
        return 0
    end
end

function _weight_probability_further(idx::Int64, reference::Int64, max_difference::Int64)::Int64
    difference = abs(idx - reference)
    if difference == max_difference
        return 16
    elseif difference == max_difference - 1
        return 13
    elseif difference == max_difference - 2
        return 10
    elseif difference == max_difference - 3
        return 7
    elseif difference == max_difference - 4
        return 4
    elseif difference == max_difference - 5
        return 1
    else
        return 0
    end
end
