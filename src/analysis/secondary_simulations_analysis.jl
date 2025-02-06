function scaling(P0_options::Vector{Float64}, number_options::Vector{Int64}, timestamp::String)
    s = readdlm("../output/$(timestamp)_scale_initial_s.csv", ',')
    s = NamedDimsArray{(:number_options, :P0)}(s)
    _plot_scaling(s, P0_options, number_options, timestamp)
end

function decision_scales(decision_steps::Vector{Float64}, time_horizons::Vector{Float64}, timestamp::String)
    s = readdlm("../output/$(timestamp)_decision_scale_s.csv", ',')
    s = NamedDimsArray{(:decision_step, :time_horizon)}(s)
    _plot_decision_scales(s, time_horizons, decision_steps, timestamp)
end
