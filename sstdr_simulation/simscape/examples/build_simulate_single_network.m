%% Build and Simulate Single Network Example
% This example demonstrates how to use the build_and_simulate_network function
% to create network configurations, build Simscape models, and run SSTDR simulations.

clear; clc; close all;

%% Setup paths
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'config'));
addpath(fullfile(parent_dir, 'functions', 'generation'));
addpath(fullfile(parent_dir, 'functions', 'network'));
addpath(fullfile(parent_dir, 'functions', 'simulation'));

fprintf('=== Build and Simulate Single Network Example ===\n');
fprintf('This example demonstrates the complete workflow from network config to simulation results.\n\n');

%% Example 1: Simple network with no faults
fprintf('--- Example 1: Simple Network (No Faults) ---\n');

% Create a simple 4-segment network with no faults
network_config_1 = create_network_config([0.1, -0.05, 0, 0.2, -0.3, 1], ...
    'dx', 1e3, ...  % 1 km per segment
    'Z0', 50, ...    % 50 ohm characteristic impedance
    'name', 'simple_no_faults');

% Create simulation configuration
sim_config_1 = create_simulation_config( ...
    'carrier_freq', 400e3, ...      
    'chip_rate', 400e3, ...         
    'fs', 3.2e6, ...     
    'pn_bits', 12, ... 
    'modulation', 'sine', ... 
    'method', 'freq', ...           % Frequency domain correlation
    'plot_results', true, ...       % Show plots
    'apply', false, ...             % Don't apply to workspace yet
    'name', 'example_1_config');

% Run the complete workflow
fprintf('Running simulation workflow...\n');
[results_1, model_info_1, sim_data_1] = build_and_simulate_network( ...
    network_config_1, sim_config_1, ...
    'save_model', true, ...        % Don't save model
    'close_model', false, ...        % Close model after simulation
    'verbose', true);

if results_1.success
    fprintf('✓ Example 1 completed successfully!\n');
    fprintf('  Simulation time: %.2f seconds\n', results_1.sim_time);
else
    fprintf('✗ Example 1 failed: %s\n', results_1.error.message);
end

% Don't do any more
return;

%% Example 2: Network with series faults
fprintf('\n--- Example 2: Network with Series Faults ---\n');

% Create a network with series faults (opens)
network_config_2 = create_network_config([0, 0.3, 0, 0.7, 0], ...
    'dx', 5.0, ...   % 5 meters per segment
    'Z0', 75, ...    % 75 ohm characteristic impedance
    'name', 'series_faults');

% Create simulation configuration with different parameters
sim_config_2 = create_simulation_config( ...
    'carrier_freq', 1e6, ...        % 1 MHz carrier
    'chip_rate', 500e3, ...         % 500 kHz chip rate
    'fs', 4e6, ...                  % 4 MHz sampling
    'method', 'both', ...           % Compare both methods
    'plot_results', true, ...       % Show plots
    'positive_only', true, ...      % Only positive time lags
    'apply', false, ...             % Don't apply to workspace yet
    'name', 'example_2_config');

% Run the complete workflow
fprintf('Running simulation workflow...\n');
[results_2, model_info_2, sim_data_2] = build_and_simulate_network( ...
    network_config_2, sim_config_2, ...
    'save_model', false, ...        % Don't save model
    'close_model', true, ...        % Close model after simulation
    'verbose', true);

if results_2.success
    fprintf('✓ Example 2 completed successfully!\n');
    fprintf('  Simulation time: %.2f seconds\n', results_2.sim_time);
    if isfield(results_2, 'peak_times') && ~isempty(results_2.peak_times)
        fprintf('  Detected peaks at: ');
        fprintf('%.3f ', results_2.peak_times * 1e6);
        fprintf('μs\n');
    end
else
    fprintf('✗ Example 2 failed: %s\n', results_2.error.message);
end

%% Example 3: Network with shunt faults
fprintf('\n--- Example 3: Network with Shunt Faults ---\n');

% Create a network with shunt faults (shorts)
network_config_3 = create_network_config([0, -0.2, 0, -0.5, 0, 0], ...
    'dx', 8.0, ...   % 8 meters per segment
    'Z0', 50, ...    % 50 ohm characteristic impedance
    'name', 'shunt_faults');

% Create simulation configuration with unmodulated signal
sim_config_3 = create_simulation_config( ...
    'modulation', 'none', ...       % Unmodulated (baseband)
    'chip_rate', 100e3, ...         % 100 kHz chip rate
    'fs', 1e6, ...                  % 1 MHz sampling
    'method', 'freq', ...           % Frequency domain correlation
    'plot_results', true, ...       % Show plots
    'apply', false, ...             % Don't apply to workspace yet
    'name', 'example_3_config');

% Run the complete workflow
fprintf('Running simulation workflow...\n');
[results_3, model_info_3, sim_data_3] = build_and_simulate_network( ...
    network_config_3, sim_config_3, ...
    'save_model', false, ...        % Don't save model
    'close_model', true, ...        % Close model after simulation
    'verbose', true);

if results_3.success
    fprintf('✓ Example 3 completed successfully!\n');
    fprintf('  Simulation time: %.2f seconds\n', results_3.sim_time);
else
    fprintf('✗ Example 3 failed: %s\n', results_3.error.message);
end

%% Example 4: Mixed faults with model saving
fprintf('\n--- Example 4: Mixed Faults (Save Model) ---\n');

% Create a network with mixed faults
network_config_4 = create_network_config([0, 0.4, -0.3, 0, 0.6, -0.1], ...
    'dx', 3.0, ...   % 3 meters per segment
    'Z0', 50, ...    % 50 ohm characteristic impedance
    'name', 'mixed_faults');

% Create simulation configuration with cosine modulation
sim_config_4 = create_simulation_config( ...
    'modulation', 'cosine', ...     % Cosine modulation
    'carrier_freq', 750e3, ...      % 750 kHz carrier
    'chip_rate', 200e3, ...         % 200 kHz chip rate
    'fs', 1.5e6, ...               % 1.5 MHz sampling
    'method', 'freq', ...           % Frequency domain correlation
    'plot_results', false, ...      % Don't show plots (for speed)
    'apply', false, ...             % Don't apply to workspace yet
    'name', 'example_4_config');

% Run the complete workflow with model saving
fprintf('Running simulation workflow with model saving...\n');
[results_4, model_info_4, sim_data_4] = build_and_simulate_network( ...
    network_config_4, sim_config_4, ...
    'save_model', true, ...         % Save model to disk
    'close_model', false, ...       % Keep model open for inspection
    'verbose', true);

if results_4.success
    fprintf('✓ Example 4 completed successfully!\n');
    fprintf('  Simulation time: %.2f seconds\n', results_4.sim_time);
    fprintf('  Model saved: %s.slx\n', results_4.model_name);
    fprintf('  Model remains open for inspection\n');
else
    fprintf('✗ Example 4 failed: %s\n', results_4.error.message);
end

%% Summary and Analysis
fprintf('\n=== Summary of All Examples ===\n');

examples = {results_1, results_2, results_3, results_4};
example_names = {'Simple (No Faults)', 'Series Faults', 'Shunt Faults', 'Mixed Faults'};

for i = 1:length(examples)
    result = examples{i};
    fprintf('%d. %s: ', i, example_names{i});
    
    if result.success
        fprintf('✓ Success (%.2f s)\n', result.sim_time);
        fprintf('   Network: %d segments, %d faults, %.1f m total\n', ...
            result.network_config.num_segments, ...
            result.network_analysis.num_faults, ...
            result.network_analysis.total_length);
        
        if isfield(result, 'peak_times') && ~isempty(result.peak_times)
            fprintf('   Peaks: %d detected\n', length(result.peak_times));
        end
    else
        fprintf('✗ Failed\n');
    end
end

%% Interactive Section
fprintf('\n=== Interactive Section ===\n');
fprintf('The examples above demonstrate different network configurations.\n');
fprintf('You can now:\n');
fprintf('1. Examine the results structures (results_1, results_2, etc.)\n');
fprintf('2. Look at the model_info structures for block details\n');
fprintf('3. Inspect the correlation results\n');
fprintf('4. Create your own network configurations\n');

% Example of creating a custom network interactively
fprintf('\n--- Create Your Own Network ---\n');
fprintf('Example: Create a custom network configuration\n');
fprintf('network_config = create_network_config([0, 0.5, -0.3, 0.2], ''dx'', 12, ''name'', ''my_network'');\n');
fprintf('sim_config = create_simulation_config(''carrier_freq'', 800e3, ''apply'', false);\n');
fprintf('[results, model_info, sim_data] = build_and_simulate_network(network_config, sim_config);\n');

fprintf('\nExample complete! Check the workspace for all results.\n'); 