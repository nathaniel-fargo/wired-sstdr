%% Test Build and Simulate Integration
% Simple test to verify the build_and_simulate_network function works correctly
% with the existing create_network_config, create_simulation_config, and other functions

clear; clc; close all;

%% Setup paths
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'config'));
addpath(fullfile(parent_dir, 'functions', 'generation'));
addpath(fullfile(parent_dir, 'functions', 'network'));
addpath(fullfile(parent_dir, 'functions', 'simulation'));

fprintf('=== Testing Build and Simulate Integration ===\n');

%% Test 1: Basic functionality test
fprintf('\n--- Test 1: Basic Functionality ---\n');

try
    % Create a simple network configuration
    network_config = create_network_config([0, 0.2, 0], ...
        'dx', 5.0, ...
        'Z0', 50, ...
        'name', 'test_network');
    
    fprintf('✓ Network configuration created successfully\n');
    
    % Create a simulation configuration
    sim_config = create_simulation_config( ...
        'carrier_freq', 500e3, ...
        'chip_rate', 100e3, ...
        'fs', 1e6, ...
        'method', 'freq', ...
        'plot_results', false, ...  % No plots for testing
        'apply', false, ...
        'name', 'test_sim_config');
    
    fprintf('✓ Simulation configuration created successfully\n');
    
    % Test the integration function
    fprintf('Running build_and_simulate_network...\n');
    [results, model_info, sim_data] = build_and_simulate_network( ...
        network_config, sim_config, ...
        'save_model', false, ...
        'close_model', true, ...
        'verbose', false);  % Reduce verbosity for testing
    
    if results.success
        fprintf('✓ Integration test passed!\n');
        fprintf('  Network: %d segments, %d faults\n', ...
            results.network_config.num_segments, ...
            results.network_analysis.num_faults);
        fprintf('  Simulation time: %.3f seconds\n', results.sim_time);
    else
        fprintf('✗ Integration test failed: %s\n', results.error.message);
    end
    
catch ME
    fprintf('✗ Test 1 failed with error: %s\n', ME.message);
    fprintf('  File: %s, Line: %d\n', ME.stack(1).file, ME.stack(1).line);
end

%% Test 2: Parameter validation
fprintf('\n--- Test 2: Parameter Validation ---\n');

try
    % Test with invalid network config (should fail gracefully)
    fprintf('Testing error handling with invalid inputs...\n');
    
    % This should fail because we're passing wrong types
    [results_bad, ~, ~] = build_and_simulate_network( ...
        'invalid_network_config', 'invalid_sim_config', ...
        'verbose', false);
    
    if ~results_bad.success
        fprintf('✓ Error handling works correctly\n');
    else
        fprintf('⚠ Error handling may not be working properly\n');
    end
    
catch ME
    % This is expected behavior
    fprintf('✓ Input validation works (caught expected error)\n');
end

%% Test 3: Configuration structure validation
fprintf('\n--- Test 3: Configuration Structure Validation ---\n');

try
    % Create configurations and check their structure
    network_config = create_network_config([0, 0.5, -0.3], 'name', 'validation_test');
    sim_config = create_simulation_config('apply', false, 'name', 'validation_test');
    
    % Check network config structure
    required_network_fields = {'name', 'type', 'num_segments', 'load_vector', ...
                              'physical', 'analysis', 'faults', 'metadata'};
    
    network_valid = true;
    for i = 1:length(required_network_fields)
        if ~isfield(network_config, required_network_fields{i})
            fprintf('✗ Missing network field: %s\n', required_network_fields{i});
            network_valid = false;
        end
    end
    
    if network_valid
        fprintf('✓ Network configuration structure is valid\n');
    end
    
    % Check simulation config structure
    required_sim_fields = {'name', 'pn_config', 'correlation_config', 'simulation_config'};
    
    sim_valid = true;
    for i = 1:length(required_sim_fields)
        if ~isfield(sim_config, required_sim_fields{i})
            fprintf('✗ Missing simulation field: %s\n', required_sim_fields{i});
            sim_valid = false;
        end
    end
    
    if sim_valid
        fprintf('✓ Simulation configuration structure is valid\n');
    end
    
catch ME
    fprintf('✗ Test 3 failed: %s\n', ME.message);
end

%% Test 4: Function availability check
fprintf('\n--- Test 4: Function Availability Check ---\n');

required_functions = {'create_network_config', 'create_simulation_config', ...
                     'build_network_model', 'run_simulation', 'build_and_simulate_network'};

all_available = true;
for i = 1:length(required_functions)
    if exist(required_functions{i}, 'file') == 2
        fprintf('✓ %s: Available\n', required_functions{i});
    else
        fprintf('✗ %s: Not found\n', required_functions{i});
        all_available = false;
    end
end

if all_available
    fprintf('✓ All required functions are available\n');
else
    fprintf('✗ Some required functions are missing\n');
end

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('Integration testing completed.\n');
fprintf('If all tests passed, the build_and_simulate_network function\n');
fprintf('is properly integrated with the existing codebase.\n');
fprintf('\nNext steps:\n');
fprintf('1. Run the full example: build_simulate_single_network.m\n');
fprintf('2. Test with your own network configurations\n');
fprintf('3. Explore the results structures for analysis\n'); 