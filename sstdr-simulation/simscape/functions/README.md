# Functions Directory Organization

The functions directory has been organized into two main categories for better code organization and maintainability:

## 📁 `network_generation/`
Contains all functions related to creating, configuring, and building network models:

### Network Configuration Functions:
- **`create_network_config.m`** - Creates network configuration structures from load vectors
- **`create_specific_network.m`** - Creates specific network patterns for testing
- **`generate_random_network.m`** - Generates random network configurations
- **`generate_network_dataset.m`** - Generates multiple random networks for datasets

### Network Model Building Functions:
- **`build_network_model.m`** - Main function to build Simscape models from network configs

### Network Utility Functions:
- **`map_load_to_resistance.m`** - Maps load values to resistance values

## 📁 `sstdr_simulation/`
Contains all functions related to SSTDR signal generation, simulation, and analysis:

### SSTDR Signal Generation:
- **`gen_pn_code.m`** - Generates PN sequences with configurable modulation
- **`configure_model_sampling.m`** - Configures Simulink model sampling parameters

### SSTDR Analysis:
- **`cross_correlate.m`** - Performs cross-correlation analysis (time/frequency domain)

## 🔧 Usage

When using these functions, make sure to add both subdirectories to your MATLAB path:

```matlab
% Add both function directories to path
addpath('functions/network_generation');
addpath('functions/sstdr_simulation');
addpath('config');
```

This path setup is automatically handled in:
- `generate_sstdr_dataset.m`
- `demo_dataset_analysis.m`
- `run_sstdr_analysis.m`
- All example scripts in `examples/`
- All test scripts in `tests/`

## 📋 Function Dependencies

### Network Generation Flow:
1. `generate_random_network()` → `create_network_config()` → network configuration
2. `build_network_model()` → Simulink model ready for simulation

### SSTDR Simulation Flow:
1. `sstdr_config()` → `gen_pn_code()` → PN sequence generation
2. `configure_model_sampling()` → model setup for SSTDR
3. `cross_correlate()` → analysis of simulation results

## 🚀 Complete Workflow Example:

```matlab
% 1. Generate network
net_config = generate_random_network('num_segments', 8);

% 2. Build model
model_name = build_network_model(net_config, 'connect_blocks', true);

% 3. Configure SSTDR
sstdr_config('default');
configure_model_sampling(model_name);

% 4. Run simulation and analysis
run_sstdr_analysis('skip', model_name);
```

This organization makes the codebase more maintainable and easier to understand by clearly separating network-related functions from SSTDR-specific functions. 