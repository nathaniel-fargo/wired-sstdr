# SSTDR Simulink/Simscape Integration Guide

This guide shows how to seamlessly integrate the enhanced SSTDR workflow with your Simscape models for comprehensive spread spectrum time domain reflectometry simulations.

## Quick Start

### 1. Basic Usage
```matlab
% Load configuration and generate PN code
sstdr_config('default');

% Configure your model's sampling parameters automatically
configure_model_sampling('my_sstdr_model');

% Run complete analysis with your Simscape model
run_sstdr_analysis('default', 'my_sstdr_model');
```

### 2. One-Command Setup
```matlab
% Load config and configure model in one step
configure_model_sampling('my_model', 'sine_fast');

% Or run complete analysis (includes model configuration)
run_sstdr_analysis('sine_fast', 'my_model');
```

### 3. Different Configurations
```matlab
% Fast sine wave modulation
run_sstdr_analysis('sine_fast', 'my_model');

% Unmodulated PN sequence
run_sstdr_analysis('unmodulated', 'my_model');

% High resolution analysis
run_sstdr_analysis('high_res', 'my_model');
```

## Simscape Model Setup

### Required Components

#### 1. Signal Source
- **Component**: `From Workspace` block
- **Variable name**: `pn_timeseries_interp` (for modulated signals) or `pn_timeseries_chip`
- **Interpolation**: Linear
- **Extrapolation**: Hold final value

#### 2. Signal Measurement
- **Component**: `To Workspace` block  
- **Variable name**: `simout` (default expected by analysis scripts)
- **Save format**: Structure with time
- **Sample time**: `sim_ts` (automatically set by configuration)
- **Maximum data points**: inf

#### 3. Transmission Line Model
Use appropriate Simscape components:
- **RF Blockset**: `TransmissionLine` blocks
- **Simscape Electrical**: Distributed parameter line models
- **Custom**: Lumped element equivalent circuits

### Model Configuration

#### Solver Settings
The configuration system automatically sets optimal solver parameters:
- **Stop Time**: Automatically calculated based on PN sequence length
- **Solver**: Optimized for the signal type (ode23t for unmodulated, ode45 for modulated)
- **Max Step Size**: Set to ensure accurate high-frequency capture

#### Sample Time Considerations
The configuration system automatically exports sampling parameters:
- **`sim_fs`**: Sampling frequency (Hz) - use for calculations
- **`sim_ts`**: Sample time (s) - set in measurement blocks  
- **`sim_decimation`**: Interpolation factor - for advanced processing

**Block Configuration:**
- Signal source sample time: Use `sim_ts` or leave as `-1` (inherited)
- Measurement block sample time: Set to `sim_ts` for consistent sampling
- Use fixed-step solvers when possible for deterministic timing

## Callback Integration

### Model Callbacks
Add these callbacks to your Simscape model for seamless operation:

#### PreLoadFcn Callback
```matlab
% Auto-load SSTDR configuration
if ~exist('sstdr_config', 'var')
    fprintf('Loading default SSTDR configuration...\n');
    sstdr_config('default');
end
```

#### PostLoadFcn Callback  
```matlab
% Verify SSTDR setup
if exist('pn_timeseries_interp', 'var')
    fprintf('✓ SSTDR signals loaded\n');
    fprintf('  Duration: %.2f ms\n', max(pn_timeseries_interp.Time)*1000);
else
    warning('⚠ SSTDR signals not found. Run sstdr_config() first.');
end
```

#### StopFcn Callback
```matlab
% Auto-run correlation analysis
if exist('out', 'var')
    fprintf('Running correlation analysis...\n');
    cross_correlate('method', 'freq');
end
```

### Block Callbacks

#### Signal Source Mask Callback
```matlab
% Update source based on current configuration
if exist('sstdr_config', 'var')
    cfg = sstdr_config;
    if strcmp(cfg.pn_config.modulation, 'none')
        set_param(gcb, 'VariableName', 'pn_timeseries_chip');
    else
        set_param(gcb, 'VariableName', 'pn_timeseries_interp');
    end
end
```

#### Measurement Block Configuration
```matlab
% Configure measurement block sample time
if exist('sim_ts', 'var')
    % Set sample time for To Workspace block
    set_param('my_model/Measurement', 'SampleTime', 'sim_ts');
    
    % Or for Scope blocks
    set_param('my_model/Scope', 'SampleTime', 'sim_ts');
    
    fprintf('✓ Measurement sample time set to %.2e s\n', evalin('base', 'sim_ts'));
else
    warning('⚠ sim_ts not found. Run sstdr_config() first.');
end
```

## Advanced Integration Patterns

### 1. Parameter Sweep Integration
```matlab
% Sweep different cable lengths
cable_lengths = [1, 5, 10, 20, 50]; % meters
results_sweep = cell(length(cable_lengths), 1);

for i = 1:length(cable_lengths)
    % Update model parameter
    set_param('my_model/Cable', 'Length', num2str(cable_lengths(i)));
    
    % Run analysis
    run_sstdr_analysis('default', 'my_model', 'save_results', true);
    
    % Store results
    results_sweep{i} = evalin('base', 'correlation_results');
end
```

### 2. Monte Carlo Analysis
```matlab
% Monte Carlo with noise variations
num_runs = 100;
noise_levels = 0.01:0.01:0.1;

for noise_level = noise_levels
    for run = 1:num_runs
        % Update noise parameters in model
        set_param('my_model/Noise', 'Variance', num2str(noise_level^2));
        
        % Run single analysis
        run_sstdr_analysis('skip', 'my_model', 'method', 'freq');
        
        % Process results...
    end
end
```

### 3. Real-time Parameter Updates
```matlab
% Function to update model parameters based on configuration
function update_simscape_model(model_name, config)
    % Update solver settings
    set_param(model_name, 'StopTime', num2str(config.simulation_config.stop_time));
    set_param(model_name, 'Solver', config.simulation_config.solver);
    set_param(model_name, 'MaxStep', num2str(config.simulation_config.max_step));
    
    % Update signal source
    if strcmp(config.pn_config.modulation, 'none')
        signal_var = 'pn_timeseries_chip';
    else
        signal_var = 'pn_timeseries_interp';  
    end
    set_param([model_name '/PN_Source'], 'VariableName', signal_var);
    
    % Update measurement settings with proper sample time
    sample_time = 1/config.pn_config.fs;
    set_param([model_name '/Measurement'], 'SampleTime', num2str(sample_time));
    set_param([model_name '/Measurement'], 'MaxDataPoints', 'inf');
    
    fprintf('✓ Model configured:\n');
    fprintf('  - Sampling: %.0f Hz (%.2e s)\n', config.pn_config.fs, sample_time);
    fprintf('  - Duration: %.3f ms\n', config.simulation_config.stop_time*1000);
    fprintf('  - Signal: %s\n', signal_var);
end
```

## Custom Analysis Functions

### Fault Location Analysis
```matlab
function fault_locations = analyze_fault_locations(results, cable_length, velocity_factor)
    % Extract reflection peaks (excluding direct signal at t=0)
    if isfield(results, 'peak_times_freq')
        reflection_times = results.peak_times_freq(results.peak_times_freq > 1e-6);
    else
        reflection_times = results.peak_times_time(results.peak_times_time > 1e-6);
    end
    
    % Calculate fault locations
    c = 3e8; % Speed of light
    v_prop = c * velocity_factor; % Propagation velocity
    fault_locations = (reflection_times * v_prop) / 2; % Two-way travel
    
    fprintf('Detected faults at:\n');
    for i = 1:length(fault_locations)
        fprintf('  %.2f m (%.1f%% of cable length)\n', ...
                fault_locations(i), 100*fault_locations(i)/cable_length);
    end
end
```

### Multi-frequency Analysis
```matlab
function run_multifrequency_analysis(model_name, frequencies)
    for i = 1:length(frequencies)
        % Generate PN code with different carrier frequency
        config = gen_pn_code('modulation', 'sine', 'carrier_freq', frequencies(i));
        
        % Run simulation
        sim_out = sim(model_name);
        assignin('base', 'out', sim_out);
        
        % Analyze results
        results = cross_correlate('method', 'freq');
        
        % Store frequency-specific results
        assignin('base', sprintf('results_%dkHz', frequencies(i)/1000), results);
    end
end
```

## Troubleshooting

### Common Issues

#### 1. "Reference code not found"
**Solution**: Ensure you run `sstdr_config()` or `gen_pn_code()` before analysis
```matlab
sstdr_config('default'); % This generates all required variables
```

#### 2. "No simulation output found"
**Solution**: Verify your measurement block saves to variable `out.simout`
```matlab
% Check if simulation output exists
if ~exist('out', 'var')
    error('Run simulation first or check measurement block configuration');
end
```

#### 3. Timing mismatches
**Solution**: Verify model sample times match configuration
```matlab
% Check sample time consistency
if exist('sim_ts', 'var')
    config_ts = evalin('base', 'sim_ts');
    model_ts = str2double(get_param('my_model/Measurement', 'SampleTime'));
    
    if abs(model_ts - config_ts) > 1e-9
        warning('Sample time mismatch: Model=%.2e, Config=%.2e', model_ts, config_ts);
        
        % Auto-fix the sample time
        set_param('my_model/Measurement', 'SampleTime', 'sim_ts');
        fprintf('✓ Fixed measurement sample time to %.2e s\n', config_ts);
    end
else
    error('Sampling parameters not found. Run sstdr_config() first.');
end
```

#### 4. Memory issues with large simulations
**Solution**: Use frequency domain correlation and limit data points
```matlab
% Use frequency domain (faster, less memory)
cross_correlate('method', 'freq', 'plot_results', false);

% Or limit measurement data points
set_param('my_model/Measurement', 'MaxDataPoints', '100000');
```

## Performance Optimization

### Best Practices
1. **Use frequency domain correlation** for long sequences (>1000 chips)
2. **Set appropriate MaxDataPoints** to avoid memory issues
3. **Use variable-step solvers** only when necessary
4. **Pre-compile models** for repeated analysis

### Benchmarking
```matlab
% Compare correlation methods
run_sstdr_analysis('default', 'my_model', 'method', 'both');

% Results will show:
% - Time domain: X ms
% - Frequency domain: Y ms  
% - Speedup: Z.Zx (freq faster)
```

## Example Workflows

### Complete Analysis Pipeline
```matlab
%% 1. Setup
sstdr_config('sine_fast');

%% 2. Simulation  
run_sstdr_analysis('sine_fast', 'sstdr_cable_model');

%% 3. Advanced Analysis
results = evalin('base', 'correlation_results');
fault_locations = analyze_fault_locations(results, 50, 0.66); % 50m cable, 0.66 VF

%% 4. Multi-configuration comparison
configs = {'default', 'sine_fast', 'unmodulated', 'high_res'};
for i = 1:length(configs)
    run_sstdr_analysis(configs{i}, 'sstdr_cable_model', 'save_results', true);
end
```

This integration approach provides a seamless workflow from model setup through advanced analysis, with automatic parameter management and flexible configuration options. 