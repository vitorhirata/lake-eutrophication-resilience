function run_sensitivity(
        P0_options::Vector{Float64}, influx::Float64, time_horizon::Float64, scenarios::Vector{Dict{Symbol, Any}}
)::String
    s = NamedDimsArray{(:P0, :type)}(zeros(length(P0_options), length(scenarios)))

    for (idx_scenario, scenario) in enumerate(scenarios)
        number_decision = compute_number_decision(time_horizon, scenario[:decision_step])
        for (idx_P0, P0) in enumerate(P0_options)
            keyword_args = filter(((k,v),) -> !(k in [:name, :decision_step]), scenario)
            s[idx_P0, idx_scenario] = _entropy(P0, influx, scenario[:decision_step], number_decision; keyword_args...)
        end
        println("Finished $(scenario[:name])")
    end

    timestamp = @sprintf("%.0f", time())
    base_filename = "../output/$(timestamp)_sensitivity_"
    writedlm("$(base_filename)s.csv",  s, ',')
    return timestamp
end
