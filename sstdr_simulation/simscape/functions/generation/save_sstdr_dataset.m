function [dataset_path, dataset_info] = save_sstdr_dataset(results, network_config, simulation_config, data_folder, verbose)
%SAVE_SSTDR_DATASET Save SSTDR simulation dataset to .mat file
%
% This function saves essential SSTDR simulation data including:
% - Network and simulation configurations
% - Raw simulation data (TX/RX signals)
% - Cross-correlation analysis results
% - Basic metadata
%
% Usage:
%   [dataset_path, dataset_info] = save_sstdr_dataset(results, network_config, simulation_config)
%   [dataset_path, dataset_info] = save_sstdr_dataset(..., data_folder, verbose)
%
% Input:
%   results           - Results structure from build_and_simulate_network
%   network_config Ω   - Network configuration structure
%   simulation_config - Simulation configuration structure
%   data_folder       - Custom data folder (optional, default: auto-generated)
%   verbose           - Display progress messages (default: true)
%
% Output:
%   dataset_path - Full path to saved dataset file
%   dataset_info - Basic information about the saved dataset

if nargin < 4
    data_folder = '';
end
if nargin < 5
    verbose = true;
end

%% Generate dataset folder and filename
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

% Generate folder name
if isempty(data_folder)
    folder_name = sprintf('sstdr_dataset_%s', timestamp);
else
    folder_name = data_folder;
end

% Create full path
current_dir = fileparts(mfilename('fullpath'));
simscape_root = fileparts(fileparts(current_dir));  % Go up from functions/generation
datasets_dir = fullfile(simscape_root, 'datasets', folder_name);

% Create directory if it doesn't exist
if ~exist(datasets_dir, 'dir')
    mkdir(datasets_dir);
    if verbose
        fprintf('Created dataset directory: %s\n', folder_name);
    end
end

% Generate filename
filename = sprintf('sstdr_data_%s.mat', timestamp);
dataset_path = fullfile(datasets_dir, filename);

%% Create simplified dataset structure
dataset = struct();

% Basic metadata
dataset.metadata = struct();
dataset.metadata.created = datestr(now);
dataset.metadata.model_name = results.model_name;
dataset.metadata.timestamp = results.timestamp;

% Network configuration
dataset.network_config = network_config;

% Simulation configuration  
dataset.simulation_config = simulation_config;

% Signal data from simulation results
if isfield(results, 'simulation_results')
    sim_results = results.simulation_results;
    
    % TX/RX signals
    if isfield(sim_results, 'tx_signal')
        dataset.tx_signal = sim_results.tx_signal;
    end
    if isfield(sim_results, 'rx_signal')
        dataset.rx_signal = sim_results.rx_signal;
    end
    
    % Cross-correlation results
    if isfield(sim_results, 'correlation_results')
        dataset.correlation_results = sim_results.correlation_results;
    end
    
    % PN code
    if isfield(sim_results, 'pn_result')
        dataset.pn_code = sim_results.pn_result;
    end
end

% Performance metrics
dataset.performance = struct();
dataset.performance.simulation_time_s = results.sim_time;

%% Save dataset
try
    save(dataset_path, 'dataset', '-v7.3');
    if verbose
        fprintf('Dataset saved to: %s\n', dataset_path);
    end
catch ME
    error('Failed to save dataset: %s', ME.message);
end

%% Create dataset info structure
dataset_info = struct();
dataset_info.path = dataset_path;
dataset_info.folder = folder_name;
dataset_info.filename = filename;
dataset_info.created = dataset.metadata.created;
dataset_info.model_name = results.model_name;

if verbose
    fprintf('Dataset saved successfully!\n');
end

end 