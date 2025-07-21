# SSTDR Simulation System

This directory contains a complete SSTDR (Spread Spectrum Time Domain Reflectometry) simulation system with enhanced frequency control and dual correlation methods.

## 🎯 Quick Start

### Simple 2-Parameter Setup
```matlab
% Just change these two parameters and run:
MODEL_NAME = 'your_model_name';  % Your Simscape model
CONFIG_CHOICE = 1;               % Configuration number (1-8)
run('examples/simple_run.m');
```

### Three-Frequency Control System

The system now uses **explicit chip rate control** with three key parameters:

1. **Chip Rate** (Hz) - Controls PN sequence timing
2. **Carrier Frequency** (Hz) - Modulation frequency  
3. **Sampling Frequency** (Hz) - Digital sampling rate

```matlab
% New frequency control pattern:
sstdr_custom_config('chip_rate', 100e3, 'carrier_freq', 100e3, 'fs', 400e3);
```

**Default Pattern**: `chip_rate = carrier_freq`, `fs = 4 × chip_rate`

## 📁 Folder Structure

```
simscape/                          # Main directory (your MATLAB working directory)
├── README.md                      # This file
├── cleanup_models.m              # Utility to clean up scattered model files
├── models/                       # Simulink models directory
│   └── sstdr_basic.slx          # Main Simscape model
├── functions/                    # Reusable functions
│   ├── generation/              # Network generation and simulation
│   │   ├── build_and_simulate_network.m # Complete workflow
│   │   ├── generate_sstdr_dataset.m     # Dataset generation
│   │   └── save_sstdr_dataset.m         # Data saving
│   ├── network/                 # Network model building
│   │   ├── build_network_model.m        # Build Simscape models
│   │   └── create_specific_network.m    # Specific network patterns
│   └── simulation/              # SSTDR simulation and analysis
│       ├── gen_pn_code.m               # PN sequence generation
│       ├── cross_correlate.m           # Dual-domain correlation
│       ├── run_simulation.m            # Simulation execution
│       └── run_simulation_analysis.m   # Complete analysis workflow
├── config/                      # Configuration files
│   ├── create_network_config.m   # Network configuration
│   ├── create_simulation_config.m # Simulation configuration
│   ├── create_dataset_config.m   # Dataset generation configuration
│   └── line_params/             # Transmission line parameters
│       └── line_test.mat        # Test transmission line data
├── examples/                    # Example scripts (start here!)
│   ├── run_single_simulation.m  # Single simulation example
│   ├── build_simulate_single_network.m # Network building example
│   └── generate_dataset_example.m      # Dataset generation example
├── datasets/                    # Generated datasets
├── tests/                       # Test scripts
└── slprj/                       # MATLAB project files
```

## 🏗️ Model Organization

All Simulink models are now organized in the `models/` subdirectory:

- **Automatic Model Management**: Models are automatically created in `models/` during simulation
- **Automatic Cleanup**: Temporary models are automatically deleted after simulation (configurable)
- **Centralized Storage**: All models are stored in one location for better organization

### Model Lifecycle:
1. **Creation**: Models are built in memory and optionally saved to `models/`
2. **Simulation**: Models run from the `models/` directory
3. **Cleanup**: Temporary models are automatically deleted after use

### Cleanup Utility:
The `cleanup_models()` function (located in `functions/network/`) can be used to clean up any scattered model files and organize them properly. Add the network functions path and run:
```matlab
addpath('functions/network');
cleanup_models();
```

## 🚀 Predefined Configurations

| Configuration | Chip Rate | Carrier | Sampling | Use Case |
|---------------|-----------|---------|----------|----------|
| `default`     | 100 kHz   | 100 kHz | 400 kHz  | Standard analysis |
| `sine_fast`   | 250 kHz   | 250 kHz | 1 MHz    | Fast processing |
| `high_res`    | 200 kHz   | 200 kHz | 800 kHz  | High resolution |
| `unmodulated` | 125 kHz   | N/A     | 500 kHz  | Baseband analysis |

```matlab
% Load predefined configuration
sstdr_config('default');        % or 'sine_fast', 'high_res', 'unmodulated'
```

## ⚙️ Custom Configurations

### Basic Custom Setup
```matlab
% Equal chip rate and carrier frequency (recommended)
sstdr_custom_config('chip_rate', 200e3, 'carrier_freq', 200e3, 'fs', 800e3);

% Unmodulated with custom chip rate
sstdr_custom_config('chip_rate', 150e3, 'modulation', 'none', 'fs', 600e3);
```

### Advanced Custom Setup
```matlab
% Different carrier and chip rates
sstdr_custom_config('chip_rate', 100e3, 'carrier_freq', 200e3, 'fs', 800e3);

% Interactive configuration
sstdr_custom_config();  % Prompts for all parameters
```

## 🔬 Key Features

### Dual Correlation Methods
- **Time Domain**: Traditional `xcorr()` approach
- **Frequency Domain**: FFT-based `IFFT(X × conj(Y))` approach
- **Performance**: Frequency domain typically 3-5× faster

### Automatic Frequency Relationships
- **Interpolation Factor**: Automatically calculated as `fs/chip_rate`
- **Solver Selection**: Optimized based on modulation type
- **Stop Time**: Calculated as `PN_length/chip_rate`

### Enhanced Signal Generation
```matlab
% The gen_pn_code function now supports:
config = gen_pn_code('chip_rate', 100e3, 'carrier_freq', 100e3, 'fs', 400e3);

% Key outputs:
% - config.chip_rate: Actual achieved chip rate
% - config.settings.interpFactor: Calculated interpolation factor
% - config.total_duration: Signal duration in seconds
```

## 📊 Frequency Relationship Guidelines

### Recommended Ratios
- **Sampling Rate**: `fs ≥ 4 × chip_rate` (minimum for good quality)
- **Carrier Frequency**: `carrier_freq ≥ chip_rate` (typical for SSTDR)
- **Nyquist Limit**: `carrier_freq < fs/2` (avoid aliasing)

### Example Configurations
```matlab
% Long Range (low frequencies)
sstdr_custom_config('chip_rate', 50e3, 'carrier_freq', 50e3, 'fs', 200e3);

% High Resolution (balanced frequencies)  
sstdr_custom_config('chip_rate', 200e3, 'carrier_freq', 200e3, 'fs', 800e3);

% High Speed (high frequencies)
sstdr_custom_config('chip_rate', 500e3, 'carrier_freq', 500e3, 'fs', 2e6);
```

## 🔧 Integration with Simscape

### Automatic Model Configuration
```matlab
% Configure your Simscape model automatically
configure_model_sampling('your_model_name');

% Run complete analysis
run_sstdr_analysis('skip', 'your_model_name');
```

### Exported Variables
The system automatically exports these variables to base workspace:
- `sim_fs`: Sampling frequency (Hz)
- `sim_ts`: Sample time (s)
- `sim_chip_rate`: Chip rate (Hz)
- `sim_decimation`: Interpolation factor
- `pn_timeseries_interp`: Interpolated PN signal for Simscape

## 📈 Performance Comparison

Run multi-configuration comparison:
```matlab
run('examples/compare_configs.m');  % Tests all predefined configurations
```

## 🧪 Testing and Validation

### Test Chip Rate Control
```matlab
run('test_chip_rate_control.m');   % Demonstrates new frequency control
```

### Test Sampling Effects  
```matlab
run('test_sampling_effect.m');     % Shows impact of sampling rate changes
```

## 🆕 What's New in Chip Rate Control

### Before (Implicit Control)
```matlab
% Old system: chip rate was fs/interpFactor
gen_pn_code('fs', 1e6, 'interpFactor', 8);  % chip_rate = 125 kHz (implicit)
```

### After (Explicit Control)
```matlab
% New system: direct chip rate specification
gen_pn_code('chip_rate', 125e3, 'fs', 500e3);  % interpFactor = 4 (calculated)
```

### Benefits
- ✅ **Explicit Control**: Direct specification of signal timing
- ✅ **Predictable**: Clear relationship between input frequencies  
- ✅ **Flexible**: Independent control of carrier and chip rates
- ✅ **Backward Compatible**: Old scripts still work
- ✅ **Automatic**: Interpolation factor calculated automatically

## 🔍 Troubleshooting

### Common Issues
1. **"Reference code not found"** → Run `sstdr_config()` first
2. **"No simulation output found"** → Check measurement block saves to `out.simout`
3. **Frequency warnings** → Adjust ratios per guidelines above


## 📚 Documentation

- **Quick Start**: `docs/QUICK_START_SCRIPTS.md`
- **Simscape Integration**: `docs/SIMSCAPE_INTEGRATION_GUIDE.md`
- **Function Reference**: See individual function headers

---

**Note**: The system maintains backward compatibility while providing enhanced control through the new chip rate parameter. All predefined configurations have been updated to use the new pattern where chip rate equals carrier frequency and sampling frequency is 4× the chip rate.

## 🆘 Troubleshooting

If you get "function not found" errors:
1. Make sure you're in the `simscape/` directory or `examples/` subdirectory
2. Use the provided example scripts (they handle paths automatically)
3. Manually add paths if needed:
   ```matlab
   addpath('functions');
   addpath('config');
   ```

## 📝 Examples of Common Tasks

### Run with Default Settings
```matlab
cd examples
simple_run  % Edit line 7 for your model
```

### Test Custom Frequency
```matlab
cd examples  
% In quick_sstdr_run.m, uncomment:
% sstdr_custom_config('carrier_freq', 500e3, 'fs', 5e6);
```

### Compare Performance
```matlab
cd examples
compare_configs  % Automatically tests 4 configurations
```

Happy analyzing! 🎉 