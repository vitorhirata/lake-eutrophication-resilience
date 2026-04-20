# Pathway Diversity in a Lake Eutrophication Model

Project implementing pathway diversity in a lake eutrophication model. The model includes comparissons of pathway
diversity with other resilience metrics, such as distance to the basin treshold and early-warning signals, and
scenarios to show pathway diversity can capture change in agency related to the policy decision-making context.

## Setup
Install julia up and download the julia version mentioned on the manifest
```
curl -fsSL https://install.julialang.org | sh
```

Open julia inside the root with:
```
julia --project=.
```

## Example of analysis
Please check the folder examples to see functions that produce results similar to the ones showed in the paper.
For example, below we show how to compute pathway diversity for multiple initial condition (and therefore the distance
to the basin threshold) and time-horizons. To run the analysis follow this steps:
```
using PathwayDiversity

# Define the model parameters
influx = 0.1 # Initial influx
decision_step = 5.0 # Decision step of 5 years
P0_options = collect(0:0.05:2.85) # Initial conditions
time_horizons = [1, 2, 4, 7] * decision_step # How long the simulation will run

# Run the simulation. This will save the result as a csv in the ./output folder
timestamp = PathwayDiversity.run_entropy(P0_options, influx, decision_step, time_horizons)

# Make the analysis and produce a plot in the ./output folder
PathwayDiversity.distance_basin_threshold(P0_options, influx, time_horizons, decision_step, timestamp)
```

See code in the examples folder for more analysis.

## Development mode

In development mode we use a folder sandbox, where you will import the library and run the project in development mode.
```
mkdir sandbox
cd sandbox
julia --project=.
```

Then run this command inside the REPL to create the Manifest:
```
] dev ../
] instantiate
```

You can then include development packages to help development.
