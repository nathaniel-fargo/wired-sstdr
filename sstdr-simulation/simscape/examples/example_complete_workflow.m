%% Complete SSTDR Dataset Workflow Example
% This script demonstrates the complete workflow:
% 1. Create configuration
% 2. Generate dataset
% 3. Analyze results

clear; clc; close all;

%% Step 1: Create Configuration
fprintf('=== Step 1: Creating Configuration ===\n');

% Create a simple configuration with default parameters
config = create_sstdr_dataset_config( ...
    'num_networks', 5, ...
    'output_dir', 'example_dataset');

%% Step 2: Generate Dataset
fprintf('\n=== Step 2: Generating Dataset ===\n');

% Generate the dataset
generate_sstdr_dataset(config);

%% Step 3: Analyze Results
fprintf('\n=== Step 3: Analyzing Results ===\n');

% Analyze with single window display
analyze_sstdr_dataset('example_dataset', 'display', 'single');

% Also create multiple windows for first 3 networks
analyze_sstdr_dataset('example_dataset', 'display', 'first_n', 'n', 3);

fprintf('\n=== Workflow Complete ===\n');
fprintf('Generated dataset in: example_dataset/\n');
fprintf('Check the figures for visualization results\n'); 