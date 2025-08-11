%% SSTDR Dataset Generation Example
% This script demonstrates how to generate datasets of random SSTDR networks
% using the dataset generation system.

clear; clc; close all;

%% Setup paths
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'config'));
addpath(fullfile(parent_dir, 'functions', 'generation'));
addpath(fullfile(parent_dir, 'functions', 'network'));
addpath(fullfile(parent_dir, 'functions', 'simulation'));

fprintf('=== SSTDR Dataset Generation Examples ===\n');
fprintf('This script demonstrates generating datasets of random SSTDR networks.\n\n');

%% Example 1: Small test dataset (5 networks)
fprintf('--- Example 1: Small Test Dataset ---\n');

% Create dataset configuration
dataset_config_1 = create_dataset_config(...
    'name', 'test_dataset', ...
    'num_networks', 5, ...
    'segments_range', [3, 6], ...
    'length_range', [20, 50], ...
    'fault_prob', 0.4, ...
    'fault_magnitude_range', [0.2, 0.6], ...
    'fault_type_prob', [0.7, 0.3], ...  % 70% open, 30% short
    'termination_range', [40, 75], ...
    'verbose', true);

% Generate dataset
fprintf('\nGenerating small test dataset...\n');
[summary_1, failed_1] = generate_sstdr_dataset(dataset_config_1);

fprintf('✓ Example 1 complete!\n');
fprintf('  Success rate: %.1f%%\n', summary_1.generation_info.success_rate * 100);
fprintf('  Total time: %.1f seconds\n', summary_1.generation_info.total_time_s);
fprintf('  Dataset folder: %s\n', summary_1.generation_info.main_folder);

%% Example 2: Medium dataset with custom simulation config
fprintf('\n--- Example 2: Medium Dataset with Custom Simulation ---\n');

% Create custom simulation config
custom_sim_config = create_simulation_config(...
    'name', 'custom_dataset_sim', ...
    'carrier_freq', 800e3, ...
    'chip_rate', 400e3, ...
    'method', 'freq', ...
    'plot_results', false, ...  % Will be forced to false anyway
    'apply', false);

% Create dataset configuration with custom simulation
dataset_config_2 = create_dataset_config(...
    'name', 'medium_dataset', ...
    'num_networks', 10, ...
    'segments_range', [4, 8], ...
    'length_range', [15, 80], ...
    'fault_prob', 0.25, ...
    'fault_magnitude_range', [0.1, 0.7], ...
    'fault_type_prob', [0.6, 0.4], ...
    'termination_range', [25, 100], ...
    'simulation_config', custom_sim_config, ...
    'verbose', true);

% Generate dataset
fprintf('\nGenerating medium dataset with custom simulation...\n');
[summary_2, failed_2] = generate_sstdr_dataset(dataset_config_2);

fprintf('✓ Example 2 complete!\n');
fprintf('  Success rate: %.1f%%\n', summary_2.generation_info.success_rate * 100);
fprintf('  Total time: %.1f seconds\n', summary_2.generation_info.total_time_s);
fprintf('  Dataset folder: %s\n', summary_2.generation_info.main_folder);

%% Example 3: Large dataset configuration (demo only - don't run)
fprintf('\n--- Example 3: Large Dataset Configuration (Demo) ---\n');

% Create configuration for a large dataset (don't actually run)
large_dataset_config = create_dataset_config(...
    'name', 'large_production_dataset', ...
    'num_networks', 1000, ...
    'segments_range', [3, 12], ...
    'length_range', [10, 120], ...
    'fault_prob', 0.3, ...
    'fault_magnitude_range', [0.05, 0.8], ...
    'fault_type_prob', [0.6, 0.4], ...
    'termination_range', [20, 120], ...
    'verbose', false, ...  % Reduce verbosity for large datasets
    'parallel', false);    % Could enable parallel processing

fprintf('Large dataset configuration created (not executed):\n');
fprintf('  Networks: %d\n', large_dataset_config.num_networks);
fprintf('  Segments: %d - %d\n', large_dataset_config.network_generation.segments_range);
fprintf('  Length: %.1f - %.1f m\n', large_dataset_config.network_generation.length_range);
fprintf('  Fault probability: %.1f%%\n', large_dataset_config.network_generation.fault_prob * 100);

% To run this large dataset, you would use:
% [summary_3, failed_3] = generate_sstdr_dataset(large_dataset_config);

%% Example 4: Resuming interrupted dataset generation
fprintf('\n--- Example 4: Resume Functionality (Demo) ---\n');

fprintf('To resume an interrupted dataset generation from network 50:\n');
fprintf('  [summary, failed] = generate_sstdr_dataset(dataset_config, ''resume_from'', 50);\n');

%% Summary of generated datasets
fprintf('\n=== Dataset Generation Summary ===\n');
fprintf('Generated datasets:\n');
fprintf('1. %s: %d/%d networks (%.1f%% success)\n', ...
    summary_1.generation_info.main_folder, ...
    summary_1.generation_info.successful_networks, ...
    summary_1.generation_info.total_networks, ...
    summary_1.generation_info.success_rate * 100);

fprintf('2. %s: %d/%d networks (%.1f%% success)\n', ...
    summary_2.generation_info.main_folder, ...
    summary_2.generation_info.successful_networks, ...
    summary_2.generation_info.total_networks, ...
    summary_2.generation_info.success_rate * 100);

%% Inspecting generated datasets
fprintf('\n--- Inspecting Generated Datasets ---\n');

% List all datasets
datasets_dir = fullfile(fileparts(pwd), 'datasets');
fprintf('All datasets in %s:\n', datasets_dir);
dataset_folders = dir(fullfile(datasets_dir, '*dataset*'));
for i = 1:length(dataset_folders)
    if dataset_folders(i).isdir
        fprintf('  %s\n', dataset_folders(i).name);
    end
end

% Show how to load a specific network's data
fprintf('\nTo load a specific network from a dataset:\n');
fprintf('  Example: Load network 001 from test_dataset\n');
fprintf('  dataset_folder = ''%s'';\n', summary_1.generation_info.main_folder);
fprintf('  network_path = fullfile(''datasets'', dataset_folder, ''network_001'', ''sstdr_data_*.mat'');\n');
fprintf('  data = load(network_path);\n');
fprintf('  network_config = data.dataset.network_config;\n');
fprintf('  load_vector = network_config.load_vector;\n');

%% Tips for large dataset generation
fprintf('\n--- Tips for Large Dataset Generation ---\n');
fprintf('For generating large datasets (>100 networks):\n');
fprintf('1. Set verbose=false to reduce output\n');
fprintf('2. Use resume_from parameter if generation is interrupted\n');
fprintf('3. Consider enabling parallel processing for faster generation\n');
fprintf('4. Monitor disk space - datasets can be large\n');
fprintf('5. Use max_retries parameter to handle occasional simulation failures\n');

fprintf('\nExample complete! Check the datasets/ folder for generated data.\n'); 