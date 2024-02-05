# Resilience metrics comparison

Project comparing resilience metrics from different models.

## Setup
Install julia up and download the latest version of julia:
```
curl -fsSL https://install.julialang.org | sh
```

Open julia inside the root with:
```
julia --project=.
```

Then run this command inside the REPL to create the Manifest:
```
] instantiate
```

Setup the sandbox, where you will import the library and run the project in development mode.
```
cd sandbox
julia --project=.
```

Then run this command inside the REPL to create the Manifest:
```
] dev ../
] instantiate
```

## Running the model
Inside the sandbox folder open the REPL with `julia --project=.`. Then, import the cookbook file with:
```
include("cookbook.jl")
```
Now you can run any function defined in the cookbook.
