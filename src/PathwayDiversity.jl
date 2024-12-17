module PathwayDiversity

using OrdinaryDiffEq, StochasticDiffEq, DiffEqNoiseProcess
using Statistics
using StatsBase
using GLM
using Loess
using Roots
using NamedDims
using Printf
using DelimitedFiles
using Plots
using StatsPlots
using Peaks
using Random
using LinearAlgebra

include("auxiliary.jl")

include("core/simulate_equation.jl")
include("core/scenarios.jl")
include("core/core.jl")
include("core/new_method.jl")

include("simulation/time_series.jl")
include("simulation/static.jl")
include("simulation/scenarios.jl")
include("simulation/internal_dynamics.jl")

include("analysis/metrics.jl")
include("analysis/analyse.jl")
include("analysis/visualisation.jl")

end
