# Lake Eutrophication Resilience

This project compares standard resilience metrics from ecology with pathway diversity. We use a lake eutrophication
model to implement this metrics.

## Setup
Install julia up and download the julia version mentioned on the manifest
```
curl -fsSL https://install.julialang.org | sh
```

Open julia inside the root with:
```
julia --project=.
```

## Example of analysis - distance to the basin threshold
In this analysis we compute pathway diversity for multiple initial condition (and therefore the distance to the basin
threshold) and time-horizons. To run the analysis follow this steps:
```
using PathwayDiversity

# Define the model parameters
influx = 0.1
decision_step = 5.0
P0_options = collect(0:0.05:2.85)
time_horizons = [1, 2, 4, 7] * decision_step # The number in the array indicates the number of decision steps

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
