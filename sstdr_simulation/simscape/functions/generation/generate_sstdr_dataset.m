function [dataset_summary, failed_networks] = generate_sstdr_dataset(dataset_config, varargin)
%GENERATE_SSTDR_DATASET Generate a complete SSTDR dataset with random networks
%
% This function generates a dataset of SSTDR simulations using randomly
% generated network configurations based on the provided dataset_config.
%
% Usage:
%   [dataset_summary, failed_networks] = generate_sstdr_dataset(dataset_config)
%   [dataset_summary, failed_networks] = generate_sstdr_dataset(dataset_config, 'param', value, ...)
%
% Input:
%   dataset_config - Configuration structure from create_dataset_config
%
% Optional Parameters:
%   resume_from    - Resume from network number (default: 1)
%   max_retries    - Maximum retries for failed networks (default: 3)
%   cleanup_models - Clean up models after each simulation (default: true)
%
% Output:
%   dataset_summary - Summary of generated dataset
%   failed_networks - List of networks that failed to simulate

%% Parse input arguments
p = inputParser;
addParameter(p, 'resume_from', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'max_retries', 3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'cleanup_models', true, @islogical);

parse(p, varargin{:});
opts = p.Results;

%% Initialize
rng('shuffle');  % Randomize seed
start_time = tic;
failed_networks = [];
successful_networks = [];

% Create main dataset folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
if isempty(dataset_config.output.data_folder)
    main_folder = sprintf('%s_%s', dataset_config.name, timestamp);
else
    main_folder = dataset_config.output.data_folder;
end

% Get simscape root directory
current_dir = fileparts(mfilename('fullpath'));
simscape_root = fileparts(fileparts(current_dir));
main_dataset_dir = fullfile(simscape_root, 'datasets', main_folder);

if ~exist(main_dataset_dir, 'dir')
    mkdir(main_dataset_dir);
end

if dataset_config.output.verbose
    fprintf('\n=== SSTDR Dataset Generation ===\n');
    fprintf('Dataset: %s\n', dataset_config.name);
    fprintf('Networks to generate: %d\n', dataset_config.num_networks);
    fprintf('Output folder: %s\n', main_folder);
    fprintf('Starting from network: %d\n', opts.resume_from);
    fprintf('\n');
end

%% Generate networks and simulate
for i = opts.resume_from:dataset_config.num_networks
    network_start_time = tic;
    network_success = false;
    retry_count = 0;
    
    if dataset_config.output.verbose
        fprintf('--- Network %d/%d ---\n', i, dataset_config.num_networks);
    end
    
    % Retry loop for failed networks
    while ~network_success && retry_count <= opts.max_retries
        try
            % Generate random network configuration
            network_config = generate_random_network(dataset_config, i);
            
            % Create simulation config for this network
            sim_config = dataset_config.base_simulation_config;
            sim_config.name = sprintf('%s_network_%03d', dataset_config.name, i);
            
            if dataset_config.output.verbose && retry_count > 0
                fprintf('  Retry %d/%d: %s\n', retry_count, opts.max_retries, network_config.name);
            elseif dataset_config.output.verbose
                fprintf('  Generating: %s\n', network_config.name);
            end
            
            % Run simulation workflow
            network_folder = sprintf('%s/network_%03d', main_folder, i);
            [results, ~, ~] = build_and_simulate_network(network_config, sim_config, ...
                'save_data', true, ...
                'data_folder', network_folder, ...
                'verbose', false, ...  % Suppress individual network verbosity
                'save_model', dataset_config.output.save_models, ...
                'close_model', opts.cleanup_models, ...
                'delete_model', opts.cleanup_models);  % Delete models when cleanup is enabled
            
            if results.success
                network_success = true;
                successful_networks(end+1) = i;
                
                network_time = toc(network_start_time);
                if dataset_config.output.verbose
                    fprintf('  ✓ Success (%.1f s): %d segments, %d faults\n', ...
                        network_time, network_config.num_segments, network_config.analysis.num_faults);
                end
            else
                retry_count = retry_count + 1;
                if dataset_config.output.verbose
                    fprintf('  ✗ Failed (attempt %d): %s\n', retry_count, results.error.message);
                end
            end
            
        catch ME
            retry_count = retry_count + 1;
            if dataset_config.output.verbose
                fprintf('  ✗ Error (attempt %d): %s\n', retry_count, ME.message);
            end
        end
    end
    
    % Record failed network
    if ~network_success
        failed_networks(end+1) = i;
        if dataset_config.output.verbose
            fprintf('  ✗ Network %d failed after %d attempts\n', i, opts.max_retries + 1);
        end
    end
    
    % Progress update
    if dataset_config.output.verbose && mod(i, 10) == 0
        elapsed_time = toc(start_time);
        avg_time_per_network = elapsed_time / (i - opts.resume_from + 1);
        remaining_networks = dataset_config.num_networks - i;
        estimated_remaining = avg_time_per_network * remaining_networks;
        
        fprintf('\n--- Progress Update ---\n');
        fprintf('Completed: %d/%d networks\n', i, dataset_config.num_networks);
        fprintf('Success rate: %.1f%%\n', length(successful_networks) / i * 100);
        fprintf('Average time per network: %.1f s\n', avg_time_per_network);
        fprintf('Estimated remaining time: %.1f minutes\n', estimated_remaining / 60);
        fprintf('\n');
    end
end

%% Create dataset summary
total_time = toc(start_time);
dataset_summary = struct();
dataset_summary.config = dataset_config;
dataset_summary.generation_info = struct();
dataset_summary.generation_info.total_time_s = total_time;
dataset_summary.generation_info.total_networks = dataset_config.num_networks;
dataset_summary.generation_info.successful_networks = length(successful_networks);
dataset_summary.generation_info.failed_networks = length(failed_networks);
dataset_summary.generation_info.success_rate = length(successful_networks) / dataset_config.num_networks;
dataset_summary.generation_info.avg_time_per_network_s = total_time / dataset_config.num_networks;
dataset_summary.generation_info.created = datestr(now);
dataset_summary.generation_info.main_folder = main_folder;
dataset_summary.generation_info.successful_network_ids = successful_networks;
dataset_summary.generation_info.failed_network_ids = failed_networks;

% Save dataset summary
summary_path = fullfile(main_dataset_dir, 'dataset_summary.mat');
save(summary_path, 'dataset_summary', '-v7.3');

if dataset_config.output.verbose
    fprintf('\n=== Dataset Generation Complete ===\n');
    fprintf('Total time: %.1f minutes\n', total_time / 60);
    fprintf('Networks generated: %d/%d (%.1f%% success)\n', ...
        length(successful_networks), dataset_config.num_networks, ...
        length(successful_networks) / dataset_config.num_networks * 100);
    fprintf('Average time per network: %.1f seconds\n', total_time / dataset_config.num_networks);
    fprintf('Failed networks: %d\n', length(failed_networks));
    if ~isempty(failed_networks)
        fprintf('Failed network IDs: %s\n', mat2str(failed_networks));
    end
    fprintf('Dataset saved to: %s\n', main_folder);
    fprintf('Summary saved to: dataset_summary.mat\n');
end

end

function network_config = generate_random_network(dataset_config, network_id)
%GENERATE_RANDOM_NETWORK Generate a random network configuration
%
% This function creates a random network configuration based on the
% parameters specified in dataset_config.

gen_params = dataset_config.network_generation;
base_config = dataset_config.base_network_config;

% Generate random network parameters
num_segments = randi(gen_params.segments_range);
total_length = gen_params.length_range(1) + ...
    rand() * (gen_params.length_range(2) - gen_params.length_range(1));

% Generate random load vector
load_vector = zeros(1, num_segments);
for i = 1:num_segments
    if rand() < gen_params.fault_prob
        % Generate fault
        fault_magnitude = gen_params.fault_magnitude_range(1) + ...
            rand() * (gen_params.fault_magnitude_range(2) - gen_params.fault_magnitude_range(1));
        
        % Determine fault type (open vs short)
        if rand() < gen_params.fault_type_prob(1)
            % Open fault (positive value)
            load_vector(i) = fault_magnitude;
        else
            % Short fault (negative value)
            load_vector(i) = -fault_magnitude;
        end
    end
    % Otherwise remains 0 (no fault)
end

% Generate random termination impedance
termination_impedance = gen_params.termination_range(1) + ...
    rand() * (gen_params.termination_range(2) - gen_params.termination_range(1));

% Create network configuration
network_name = sprintf('random_network_%03d', network_id);
network_config = create_network_config(load_vector, ...
    'name', network_name, ...
    'dx', total_length / num_segments, ...
    'termination_impedance', termination_impedance);

% Override with base config parameters if provided
if isfield(base_config, 'physical')
    if isfield(base_config.physical, 'Z0')
        network_config.physical.Z0 = base_config.physical.Z0;
    end
    if isfield(base_config.physical, 'velocity')
        network_config.physical.velocity = base_config.physical.velocity;
    end
end

end 