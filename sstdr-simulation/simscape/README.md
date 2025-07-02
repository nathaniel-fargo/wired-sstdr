# SSTDR Simscape Simulation - Organized Structure

This folder contains a complete SSTDR (Spread Spectrum Time Domain Reflectometry) simulation and analysis workflow for MATLAB/Simscape.

## 📁 Folder Organization

```
simscape/                          # Main directory (your MATLAB working directory)
├── README.md                      # This file
├── sstdr_basic.slx               # Main Simscape model
├── run_sstdr_analysis.m          # Main workflow script
├── line_test.mat                 # Test data
├── functions/                    # Reusable functions
│   ├── gen_pn_code.m            # PN sequence generation
│   ├── cross_correlate.m        # Correlation analysis
│   ├── configure_model_sampling.m # Model configuration
│   └── line_params.m            # Line parameter utilities
├── config/                      # Configuration files
│   ├── sstdr_config.m           # Predefined configurations
│   └── sstdr_custom_config.m    # Custom configuration builder
├── examples/                    # Example scripts (start here!)
│   ├── simple_run.m             # Minimal 2-parameter setup
│   ├── quick_sstdr_run.m        # Full-featured script
│   ├── compare_configs.m        # Multi-configuration comparison
│   └── cross_corr_pure.m        # Pure MATLAB correlation demo
├── docs/                        # Documentation
│   ├── QUICK_START_SCRIPTS.md   # How to use example scripts
│   └── SIMSCAPE_INTEGRATION_GUIDE.md # Integration guide
├── archive/                     # Backup files
└── slprj/                       # MATLAB project files
```

## 🚀 Quick Start

### Option 1: Super Simple (2 parameters)
```matlab
cd examples
open simple_run.m
% Edit MODEL_NAME and CONFIG_CHOICE
% Hit run!
```

### Option 2: Full Featured
```matlab
cd examples
open quick_sstdr_run.m
% Edit MODEL_NAME and uncomment one configuration
% Hit run!
```

### Option 3: Compare Multiple Configurations
```matlab
cd examples
open compare_configs.m
% Edit MODEL_NAME (optional)
% Hit run!
```

## 📋 Available Configurations

1. **Default** - 100 kHz carrier, 1 MHz sampling
2. **Fast** - 250 kHz carrier, 2 MHz sampling  
3. **High Resolution** - 100 kHz carrier, 5 MHz sampling
4. **Unmodulated** - No carrier, 1 MHz sampling
5. **Custom** - Build your own with `sstdr_custom_config()`

## 🔧 Advanced Usage

### From Root Directory (simscape/)
```matlab
% Complete workflow
run_sstdr_analysis('default', 'my_model');

% Just correlation analysis
cross_correlate('method', 'both');

% Configure model
configure_model_sampling('my_model');
```

### From Examples Directory
All example scripts automatically set up the necessary paths.

## 📖 Documentation

- **Quick Start Guide**: `docs/QUICK_START_SCRIPTS.md`
- **Integration Guide**: `docs/SIMSCAPE_INTEGRATION_GUIDE.md`

## 🏗️ Path Management

Scripts automatically add the necessary directories to your MATLAB path:
- `functions/` - Core functionality
- `config/` - Configuration management

No manual path setup required when using the provided scripts!

## 🎯 Key Features

- **Dual Correlation Methods**: Time domain (xcorr) and frequency domain (FFT)
- **Flexible Modulation**: Sine, cosine, or unmodulated PN sequences
- **Auto Configuration**: Predefined setups for common use cases
- **Performance Comparison**: Speed benchmarking between methods
- **Simscape Integration**: Automatic model configuration
- **Complete Documentation**: Guides for every use case

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