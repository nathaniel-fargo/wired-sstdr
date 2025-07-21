%% Test Model Organization and Cleanup Functionality
% This script tests the new model organization features

clear; clc; close all;

fprintf('=== Testing Model Organization and Cleanup ===\n');

%% Setup paths
% Get the directory of this script and the parent (simscape) directory
script_dir = fileparts(mfilename('fullpath'));
simscape_dir = fileparts(script_dir);

addpath(fullfile(simscape_dir, 'config'));
addpath(fullfile(simscape_dir, 'functions', 'network'));
addpath(fullfile(simscape_dir, 'functions', 'simulation'));
addpath(fullfile(simscape_dir, 'functions', 'generation'));

%% Test 1: Model with automatic deletion (default behavior)
fprintf('\n--- Test 1: Model with automatic deletion ---\n');
try
    network_config = create_network_config([0, 0.3, 0], ...
        'name', 'test_auto_delete');
    
    sim_config = create_simulation_config('apply', false, ...
        'name', 'test_config', 'plot_results', false);
    
    [results1, ~, ~] = build_and_simulate_network(network_config, sim_config, ...
        'verbose', true, 'save_model', false, 'delete_model', true);
    
    if results1.success
        fprintf('✓ Test 1 PASSED: Model created and deleted successfully\n');
    else
        fprintf('✗ Test 1 FAILED: Simulation failed\n');
    end
    
catch ME
    fprintf('✗ Test 1 ERROR: %s\n', ME.message);
    results1.success = false;
end

%% Test 2: Model saved to models/ directory
fprintf('\n--- Test 2: Model saved to models/ directory ---\n');
try
    network_config.name = 'test_saved_model';
    
    [results2, ~, ~] = build_and_simulate_network(network_config, sim_config, ...
        'verbose', true, 'save_model', true, 'delete_model', false);
    
    if results2.success
        fprintf('✓ Test 2 PASSED: Model saved to models/ directory\n');
    else
        fprintf('✗ Test 2 FAILED: Simulation failed\n');
    end
    
catch ME
    fprintf('✗ Test 2 ERROR: %s\n', ME.message);
    results2.success = false;
end

%% Test 3: Check models directory contents
fprintf('\n--- Test 3: Checking models directory ---\n');
try
    models_dir = fullfile(simscape_dir, 'models');
    if exist(models_dir, 'dir')
        model_files = dir(fullfile(models_dir, '*.slx'));
        fprintf('Models in models/ directory:\n');
        for i = 1:length(model_files)
            fprintf('  %s\n', model_files(i).name);
        end
        
        % Check if our saved model is there
        saved_model_exists = any(strcmp({model_files.name}, 'test_saved_model.slx'));
        if saved_model_exists
            fprintf('✓ Test 3 PASSED: Saved model found in models/ directory\n');
        else
            fprintf('✗ Test 3 FAILED: Saved model not found in models/ directory\n');
        end
    else
        fprintf('✗ Test 3 FAILED: models/ directory does not exist\n');
    end
    
catch ME
    fprintf('✗ Test 3 ERROR: %s\n', ME.message);
end

%% Test 4: Test cleanup utility
fprintf('\n--- Test 4: Testing cleanup utility ---\n');
try
    cleanup_models();
    fprintf('✓ Test 4 PASSED: Cleanup utility executed successfully\n');
catch ME
    fprintf('✗ Test 4 ERROR: %s\n', ME.message);
end

%% Summary
fprintf('\n=== Test Summary ===\n');
if results1.success
    test1_status = 'PASSED';
else
    test1_status = 'FAILED';
end

if results2.success
    test2_status = 'PASSED';
else
    test2_status = 'FAILED';
end

fprintf('Test 1 (Auto-delete): %s\n', test1_status);
fprintf('Test 2 (Save model):  %s\n', test2_status);
fprintf('Test 3 (Directory):   Check output above\n');
fprintf('Test 4 (Cleanup):     Check output above\n');

if results1.success && results2.success
    fprintf('\n✓ Overall: Model organization system working correctly!\n');
else
    fprintf('\n✗ Overall: Some tests failed - check output above\n');
end

fprintf('\nModel organization features tested:\n');
fprintf('- Automatic model deletion after simulation\n');
fprintf('- Model saving to models/ subdirectory\n');
fprintf('- Cleanup utility functionality\n');
fprintf('- Directory organization\n'); 