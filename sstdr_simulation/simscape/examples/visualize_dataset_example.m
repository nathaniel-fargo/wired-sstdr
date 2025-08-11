%% SSTDR Dataset Visualization Examples
% This script demonstrates how to use the visualize_sstdr_dataset function
% to analyze and visualize generated SSTDR datasets

clear; clc; close all;

%% Setup paths
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'functions', 'generation'));
addpath(fullfile(parent_dir, 'config'));

fprintf('=== SSTDR Dataset Visualization Examples ===\n\n');

%% Find available datasets
datasets_dir = fullfile(parent_dir, 'datasets');
if ~isfolder(datasets_dir)
    error('No datasets folder found. Please generate some datasets first using generate_dataset_example.m');
end

% Get list of dataset folders
dataset_list = dir(datasets_dir);
dataset_list = dataset_list([dataset_list.isdir] & ~ismember({dataset_list.name}, {'.', '..'}));

if isempty(dataset_list)
    error('No datasets found. Please generate some datasets first using generate_dataset_example.m');
end

fprintf('Available datasets:\n');
for i = 1:length(dataset_list)
    fprintf('  %d. %s\n', i, dataset_list(i).name);
end

% Use a dataset with actual networks
selected_dataset = '';
dataset_path = '';

% Find a dataset with networks
for i = length(dataset_list):-1:1  % Start from most recent
    test_path = fullfile(datasets_dir, dataset_list(i).name);
    network_dirs = dir(fullfile(test_path, 'network_*'));
    network_dirs = network_dirs([network_dirs.isdir]);
    
    if ~isempty(network_dirs)
        selected_dataset = dataset_list(i).name;
        dataset_path = test_path;
        break;
    end
end

if isempty(selected_dataset)
    error('No datasets with networks found. Please generate datasets first using generate_dataset_example.m');
end

fprintf('\nUsing dataset: %s\n\n', selected_dataset);

%% Example 1: Basic Dataset Overview
fprintf('--- Example 1: Dataset Overview and Statistics ---\n');

try
    [stats, fig_handles] = visualize_sstdr_dataset(dataset_path, ...
        'verbose', true, ...
        'show_summary', true);
    
    fprintf('✓ Dataset overview complete\n');
    fprintf('  Generated %d figure(s)\n', length(fig_handles));
    fprintf('  Dataset contains %d networks\n', stats.num_networks);
    
catch ME
    fprintf('✗ Example 1 failed: %s\n', ME.message);
end

%% Example 2: Visualize Specific Network
fprintf('\n--- Example 2: Specific Network Visualization ---\n');

try
    % Visualize network 1 with detailed analysis
    [~, fig_handles2] = visualize_sstdr_dataset(dataset_path, ...
        'network_id', 1, ...
        'verbose', true, ...
        'show_summary', false);
    
    fprintf('✓ Network 1 visualization complete\n');
    fprintf('  Generated %d figure(s)\n', length(fig_handles2));
    
catch ME
    fprintf('✗ Example 2 failed: %s\n', ME.message);
end

%% Example 3: Save Plots to Files
fprintf('\n--- Example 3: Save Visualization Plots ---\n');

try
    % Create a temporary visualization with saved plots
    [~, fig_handles3] = visualize_sstdr_dataset(dataset_path, ...
        'network_id', 1, ...
        'save_plots', true, ...
        'plot_format', 'png', ...
        'verbose', false);
    
    fprintf('✓ Plots saved to network folder\n');
    fprintf('  Generated and saved %d figure(s)\n', length(fig_handles3));
    
catch ME
    fprintf('✗ Example 3 failed: %s\n', ME.message);
end

%% Example 4: Visualize Multiple Networks (if dataset is small)
if exist('stats', 'var') && isfield(stats, 'num_networks') && stats.num_networks <= 3
    fprintf('\n--- Example 4: Multiple Network Visualization ---\n');
    
    try
        [~, fig_handles4] = visualize_sstdr_dataset(dataset_path, ...
            'show_all', true, ...
            'verbose', false);
        
        fprintf('✓ All networks visualization complete\n');
        fprintf('  Generated %d figure(s) for %d networks\n', length(fig_handles4), stats.num_networks);
        
    catch ME
        fprintf('✗ Example 4 failed: %s\n', ME.message);
    end
else
    fprintf('\n--- Example 4: Skipped (dataset too large or stats unavailable) ---\n');
    if exist('stats', 'var') && isfield(stats, 'num_networks')
        fprintf('Dataset has %d networks. Skipping show_all visualization.\n', stats.num_networks);
    else
        fprintf('Dataset statistics not available.\n');
    end
end

%% Example 5: Direct File Visualization
fprintf('\n--- Example 5: Direct Network File Visualization ---\n');

try
    % Find the first network file
    network_files = dir(fullfile(dataset_path, 'network_*', '*.mat'));
    if ~isempty(network_files)
        network_file_path = fullfile(network_files(1).folder, network_files(1).name);
        
        [~, fig_handles5] = visualize_sstdr_dataset(network_file_path, ...
            'verbose', true, ...
            'show_summary', false);
        
        fprintf('✓ Direct file visualization complete\n');
        fprintf('  Generated %d figure(s)\n', length(fig_handles5));
    else
        fprintf('✗ No network files found in dataset\n');
    end
    
catch ME
    fprintf('✗ Example 5 failed: %s\n', ME.message);
end

%% Example 6: Demonstrate Statistics Analysis
fprintf('\n--- Example 6: Detailed Statistics Analysis ---\n');

try
    % Get detailed statistics without creating plots
    [detailed_stats, ~] = visualize_sstdr_dataset(dataset_path, ...
        'verbose', false, ...
        'show_summary', false);
    
    % Display custom analysis
    fprintf('=== Custom Dataset Analysis ===\n');
    fprintf('Total networks analyzed: %d\n', detailed_stats.num_networks);
    
    if detailed_stats.num_networks > 0
        fprintf('\nNetwork Complexity:\n');
        fprintf('  Average segments: %.1f (±%.1f)\n', detailed_stats.segments.mean, detailed_stats.segments.std);
        fprintf('  Segment range: %d - %d\n', detailed_stats.segments.min, detailed_stats.segments.max);
        
        fprintf('\nPhysical Characteristics:\n');
        fprintf('  Average length: %.1f m (±%.1f m)\n', detailed_stats.length.mean, detailed_stats.length.std);
        fprintf('  Length range: %.1f - %.1f m\n', detailed_stats.length.min, detailed_stats.length.max);
        
        fprintf('\nFault Characteristics:\n');
        fprintf('  Total faults in dataset: %d\n', detailed_stats.faults.total);
        fprintf('  Average faults per network: %.1f (±%.1f)\n', ...
            detailed_stats.faults.per_network_mean, detailed_stats.faults.per_network_std);
        fprintf('  Fault distribution: %d series, %d shunt\n', ...
            detailed_stats.fault_types.series, detailed_stats.fault_types.shunt);
        
        if detailed_stats.fault_types.series + detailed_stats.fault_types.shunt > 0
            series_percent = detailed_stats.fault_types.series / ...
                (detailed_stats.fault_types.series + detailed_stats.fault_types.shunt) * 100;
            fprintf('  Series fault percentage: %.1f%%\n', series_percent);
        end
        
        fprintf('\nTermination Characteristics:\n');
        fprintf('  Average termination: %.1f Ω (±%.1f Ω)\n', ...
            detailed_stats.termination.mean, detailed_stats.termination.std);
        fprintf('  Termination range: %.1f - %.1f Ω\n', ...
            detailed_stats.termination.min, detailed_stats.termination.max);
    end
    
    fprintf('✓ Detailed statistics analysis complete\n');
    
catch ME
    fprintf('✗ Example 6 failed: %s\n', ME.message);
end

%% Summary
fprintf('\n=== Visualization Examples Summary ===\n');
total_figures = 0;
if exist('fig_handles', 'var'), total_figures = total_figures + length(fig_handles); end
if exist('fig_handles2', 'var'), total_figures = total_figures + length(fig_handles2); end
if exist('fig_handles3', 'var'), total_figures = total_figures + length(fig_handles3); end
if exist('fig_handles4', 'var'), total_figures = total_figures + length(fig_handles4); end
if exist('fig_handles5', 'var'), total_figures = total_figures + length(fig_handles5); end

fprintf('Dataset visualized: %s\n', selected_dataset);
fprintf('Total figures created: %d\n', total_figures);
fprintf('\nVisualization Features Demonstrated:\n');
fprintf('✓ Dataset overview and statistics\n');
fprintf('✓ Network topology diagrams\n');
fprintf('✓ Signal plots (TX/RX)\n');
fprintf('✓ Cross-correlation analysis\n');
fprintf('✓ Fault analysis and visualization\n');
fprintf('✓ Plot saving functionality\n');
fprintf('✓ Multiple visualization modes\n');

fprintf('\n--- Tips for Using the Visualization Tool ---\n');
fprintf('1. Use visualize_sstdr_dataset(dataset_folder) for quick overview\n');
fprintf('2. Add ''network_id'', N to focus on specific network\n');
fprintf('3. Use ''save_plots'', true to save visualizations\n');
fprintf('4. Set ''verbose'', false to reduce output for large datasets\n');
fprintf('5. Use ''show_all'', true for small datasets (<10 networks)\n');
fprintf('6. Direct file paths work: visualize_sstdr_dataset(''path/to/network.mat'')\n');

fprintf('\nVisualization examples complete!\n'); 