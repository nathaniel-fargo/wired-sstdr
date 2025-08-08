### SSTDR Simulation Workspace

Goal: generate labeled SSTDR data for learning on branching networks. This folder contains two tracks:

- Literature reproductions (`edun-et-al`, `hunter-et-al`, `evan`) kept mostly intact.
- A Simscape-based, programmatic data generator under `simscape/` that builds cable networks, runs simulations, and saves datasets end-to-end.

### Simscape pipeline

1) Configure

```
cfgNet = simscape.config.create_network_config();
cfgSim = simscape.config.create_simulation_config();
cfgData = simscape.config.create_dataset_config();
```

2) Build and simulate one network

```
mdl = simscape.functions.network.build_network_model(cfgNet);
simOut = simscape.functions.simulation.run_simulation(mdl, cfgSim);
simscape.functions.simulation.run_simulation_analysis(simOut);
```

3) Generate a dataset programmatically

```
ds = simscape.functions.generation.generate_sstdr_dataset(cfgNet, cfgSim, cfgData);
simscape.functions.generation.save_sstdr_dataset(ds, 'simscape/datasets');
```

See `simscape/examples/` for end-to-end scripts (`generate_dataset_example.m`, `build_simulate_single_network.m`). Tests under `simscape/tests/` exercise the build and parameterization.

Notes:
- Models are saved under `simscape/models/` and cached in `slprj/` by MATLAB.
- `machine_learning/spec.md` outlines the intended dataset structure for training.
- The pipeline is functional but not fully validated against hardware; treat as a starting point.


