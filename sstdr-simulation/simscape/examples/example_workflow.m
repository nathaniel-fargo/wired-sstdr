%% SSTDR Dataset Generation Example Workflow
% This example demonstrates the complete workflow with the updated folder structure

clear; clc; close all;

% Add paths to dataset generation functions
addpath('../functions/dataset_generation');

%% Example 1: Default Configuration
fprintf('=== Example 1: Default Configuration ===\n');

% Create default configuration
% - Datasets automatically go to datasets/ folder with timestamp
% - Line parameters automatically loaded from config/line_params/
config1 = create_sstdr_dataset_config( ...
    'num_networks', 5);

% Generate dataset
generate_sstdr_dataset(config1);

% Analyze results
analyze_sstdr_dataset(config1.dataset.output_dir, 'display', 'single');

%% Example 2: Custom Configuration with Named Output
fprintf('\n=== Example 2: Custom Configuration ===\n');

% Create custom configuration with specific output directory
config2 = create_sstdr_dataset_config( ...
    'chip_rate', 10e6, ...
    'fs', 40e6, ...
    'pn_bits', 12, ...
    'num_networks', 3, ...
    'output_dir', 'my_custom_dataset');  % Will be placed in datasets/my_custom_dataset

% Generate dataset
generate_sstdr_dataset(config2);

% Analyze with multiple windows
analyze_sstdr_dataset(config2.dataset.output_dir, 'display', 'multiple');

%% Example 3: Different Line Parameters
fprintf('\n=== Example 3: Different Line Parameters ===\n');

% If you have other line parameter files in config/line_params/
% Just specify the filename - path is automatically handled
config3 = create_sstdr_dataset_config( ...
    'transmission_line_file', 'line_test.mat', ...  % Automatically becomes config/line_params/line_test.mat
    'num_networks', 2, ...
    'output_dir', 'line_test_dataset');

% Generate dataset
generate_sstdr_dataset(config3);

% Analyze first network only
analyze_sstdr_dataset(config3.dataset.output_dir, 'display', 'first_n', 'n', 1);

%% Show final structure
fprintf('\n=== Final Directory Structure ===\n');
fprintf('Datasets created in:\n');
fprintf('  %s\n', config1.dataset.output_dir);
fprintf('  %s\n', config2.dataset.output_dir);
fprintf('  %s\n', config3.dataset.output_dir);
fprintf('\nLine parameters loaded from:\n');
fprintf('  %s\n', config1.network.transmission_line_file);

fprintf('\nWorkflow complete!\n'); 