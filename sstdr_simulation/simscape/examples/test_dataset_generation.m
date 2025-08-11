%% Test Dataset Generation System
% This script tests the complete dataset generation workflow

clear; clc; close all;

%% Setup paths
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'config'));
addpath(fullfile(parent_dir, 'functions', 'generation'));
addpath(fullfile(parent_dir, 'functions', 'network'));
addpath(fullfile(parent_dir, 'functions', 'simulation'));

fprintf('=== Testing Dataset Generation System ===\n\n');

%% Test 1: Basic dataset generation
fprintf('--- Test 1: Basic Dataset Generation ---\n');

dataset_config = create_dataset_config(...
    'name', 'test_basic', ...
    'num_networks', 3, ...
    'segments_range', [3, 5], ...
    'length_range', [20, 40], ...
    'fault_prob', 0.4, ...
    'fault_magnitude_range', [0.2, 0.6], ...
    'fault_type_prob', [0.6, 0.4], ...
    'termination_range', [40, 80], ...
    'verbose', false);  % Reduce verbosity for test

[summary, failed] = generate_sstdr_dataset(dataset_config);

% Verify results
assert(summary.generation_info.total_networks == 3, 'Should generate 3 networks');
assert(summary.generation_info.success_rate == 1.0, 'Should have 100% success rate');
assert(isempty(failed), 'Should have no failed networks');

fprintf('✓ Basic dataset generation test passed\n');
fprintf('  Generated: %d/%d networks\n', summary.generation_info.successful_networks, summary.generation_info.total_networks);
fprintf('  Success rate: %.1f%%\n', summary.generation_info.success_rate * 100);

%% Test 2: Verify dataset structure
fprintf('\n--- Test 2: Dataset Structure Verification ---\n');

dataset_folder = summary.generation_info.main_folder;
datasets_dir = fullfile(fileparts(pwd), 'datasets', dataset_folder);

% Check main dataset folder exists
assert(exist(datasets_dir, 'dir') == 7, 'Main dataset folder should exist');

% Check dataset summary exists
summary_file = fullfile(datasets_dir, 'dataset_summary.mat');
assert(exist(summary_file, 'file') == 2, 'Dataset summary file should exist');

% Check individual network folders
for i = 1:3
    network_folder = fullfile(datasets_dir, sprintf('network_%03d', i));
    assert(exist(network_folder, 'dir') == 7, sprintf('Network %d folder should exist', i));
    
    % Check network data file exists
    network_files = dir(fullfile(network_folder, 'sstdr_data_*.mat'));
    assert(length(network_files) == 1, sprintf('Network %d should have exactly one data file', i));
end

fprintf('✓ Dataset structure verification passed\n');
fprintf('  Main folder: %s\n', dataset_folder);
fprintf('  Summary file: ✓\n');
fprintf('  Network folders: ✓ (3/3)\n');

%% Test 3: Verify data content
fprintf('\n--- Test 3: Data Content Verification ---\n');

% Load first network data
network_file = fullfile(datasets_dir, 'network_001', 'sstdr_data_*.mat');
network_files = dir(network_file);
network_data = load(fullfile(network_files(1).folder, network_files(1).name));

% Check essential fields exist
required_fields = {'metadata', 'network_config', 'simulation_config', 'tx_signal', 'rx_signal', 'correlation_results', 'performance'};
dataset_fields = fieldnames(network_data.dataset);

for i = 1:length(required_fields)
    assert(ismember(required_fields{i}, dataset_fields), sprintf('Field %s should exist', required_fields{i}));
end

% Check network config has load vector and termination impedance
assert(isfield(network_data.dataset.network_config, 'load_vector'), 'Network config should have load_vector');
assert(isfield(network_data.dataset.network_config.physical, 'termination_impedance'), 'Network config should have termination_impedance');

% Check load vector is valid
load_vector = network_data.dataset.network_config.load_vector;
assert(all(load_vector >= -1 & load_vector <= 1), 'Load vector values should be in [-1, 1]');

% Check termination impedance is reasonable
term_impedance = network_data.dataset.network_config.physical.termination_impedance;
assert(term_impedance >= 40 && term_impedance <= 80, 'Termination impedance should be in expected range');

fprintf('✓ Data content verification passed\n');
fprintf('  Required fields: ✓ (%d/%d)\n', length(required_fields), length(required_fields));
fprintf('  Load vector: ✓ (length %d)\n', length(load_vector));
fprintf('  Termination impedance: ✓ (%.1f Ω)\n', term_impedance);

%% Test 4: Custom simulation config
fprintf('\n--- Test 4: Custom Simulation Config ---\n');

custom_sim_config = create_simulation_config(...
    'name', 'test_custom', ...
    'carrier_freq', 500e3, ...
    'chip_rate', 250e3, ...
    'method', 'freq', ...
    'plot_results', false, ...
    'apply', false);

dataset_config_custom = create_dataset_config(...
    'name', 'test_custom_sim', ...
    'num_networks', 2, ...
    'segments_range', [4, 6], ...
    'simulation_config', custom_sim_config, ...
    'verbose', false);

[summary_custom, failed_custom] = generate_sstdr_dataset(dataset_config_custom);

assert(summary_custom.generation_info.success_rate == 1.0, 'Custom simulation should succeed');
assert(isempty(failed_custom), 'Should have no failed networks with custom config');

fprintf('✓ Custom simulation config test passed\n');
fprintf('  Generated: %d/%d networks\n', summary_custom.generation_info.successful_networks, summary_custom.generation_info.total_networks);

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('All tests passed successfully!\n');
fprintf('✓ Basic dataset generation\n');
fprintf('✓ Dataset structure verification\n');
fprintf('✓ Data content verification\n');
fprintf('✓ Custom simulation config\n');
fprintf('\nDataset generation system is working correctly.\n');

%% Cleanup info
fprintf('\n--- Generated Test Datasets ---\n');
fprintf('1. %s (3 networks)\n', summary.generation_info.main_folder);
fprintf('2. %s (2 networks)\n', summary_custom.generation_info.main_folder);
fprintf('These can be found in the datasets/ folder.\n'); 