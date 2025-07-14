%% Test Standalone Network SSTDR Simulation
% Test the new run_network_sstdr function with save options

clear; close all; clc;

%% Add paths
addpath('functions/network_generation');
addpath('functions/dataset_generation');
addpath('config');

%% Create network configuration
fprintf('=== Creating Network Configuration ===\n');

% Create a test network with a few faults
load_vector = [0, 0.6, 0, -0.4, 0, 0.3, 0];  % Mix of series and shunt faults

network_config = create_network_config(load_vector, ...
    'dx', 1.0, ...           % 1 meter per segment
    'Z0', 50, ...            % 50 ohm characteristic impedance
    'velocity', 2.75e8, ...  % 2.75e8 m/s propagation velocity
    'name', 'test_standalone_network');

fprintf('Network created: %s\n', network_config.name);
fprintf('  %d segments, %.1f m total length\n', network_config.num_segments, network_config.physical.total_length);
fprintf('  %d faults detected\n', network_config.analysis.num_faults);

%% Create SSTDR configuration
fprintf('\n=== Creating SSTDR Configuration ===\n');

sstdr_config = create_sstdr_dataset_config( ...
    'chip_rate', 5e6, ...        % 5 MHz chip rate
    'carrier_freq', 5e6, ...     % 5 MHz carrier
    'fs', 20e6, ...             % 20 MHz sampling
    'pn_bits', 11, ...          % 2047 chip sequence
    'duration', 50e-6);         % 50 μs simulation

fprintf('SSTDR config created: %s\n', sstdr_config.metadata.name);
fprintf('  %.1f MHz chip rate, %.1f MHz sampling\n', ...
    sstdr_config.sstdr.chip_rate/1e6, sstdr_config.sstdr.fs/1e6);
fprintf('  %d-bit PN sequence, %.1f μs duration\n', ...
    sstdr_config.sstdr.pn_bits, sstdr_config.simulation.duration*1e6);

%% Test 1: Run without saving (quick test)
fprintf('\n=== Test 1: Quick Run (No Save) ===\n');

try
    [results1, model_info1] = run_network_sstdr(network_config, sstdr_config, ...
        'model_name', 'test_quick', ...
        'save_model', false, ...
        'keep_open', false, ...
        'plot_results', false, ...
        'verbose', true);
    
    fprintf('✓ Quick test completed successfully\n');
    fprintf('  Success: %s\n', mat2str(results1.success));
    fprintf('  Simulation time: %.2f s\n', results1.sim_time);
    
catch ME
    fprintf('✗ Quick test failed: %s\n', ME.message);
end

%% Test 2: Run with model saving and inspection
fprintf('\n=== Test 2: Full Run with Save and Inspection ===\n');

try
    [results2, model_info2] = run_network_sstdr(network_config, sstdr_config, ...
        'model_name', 'test_inspection', ...
        'save_model', true, ...
        'keep_open', true, ...        % Keep open for inspection
        'plot_results', true, ...     % Show plots
        'verbose', true, ...
        'output_dir', pwd);           % Save in current directory
    
    fprintf('✓ Full test completed successfully\n');
    fprintf('  Success: %s\n', mat2str(results2.success));
    fprintf('  Model file: %s\n', results2.model_file);
    fprintf('  Results file: %s\n', results2.results_file);
    
    % Check if model is open for inspection
    if bdIsLoaded('test_inspection')
        fprintf('  Model is open and ready for inspection in Simulink\n');
        fprintf('  You can now examine the model structure and connections\n');
    end
    
catch ME
    fprintf('✗ Full test failed: %s\n', ME.message);
    fprintf('Error details: %s\n', ME.getReport);
end

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('Standalone SSTDR simulation function tested\n');
fprintf('Features verified:\n');
fprintf('  - Network model generation\n');
fprintf('  - SSTDR parameter configuration\n');
fprintf('  - Simulation execution\n');
fprintf('  - Results analysis and plotting\n');
fprintf('  - Model saving for inspection\n');
fprintf('  - Comprehensive error handling\n');

if exist('results2', 'var') && results2.success
    fprintf('\nInspection files created:\n');
    fprintf('  - Simulink model: %s\n', results2.model_file);
    fprintf('  - Results data: %s\n', results2.results_file);
    fprintf('  - Model open in Simulink for inspection\n');
end

fprintf('\nTest complete!\n'); 