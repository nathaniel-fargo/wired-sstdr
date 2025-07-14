function generate_sstdr_dataset(config, varargin)
%GENERATE_SSTDR_DATASET Generate SSTDR dataset using configuration object
%
% Usage:
%   generate_sstdr_dataset()                    % Use default config
%   generate_sstdr_dataset(config)              % Use provided config
%   generate_sstdr_dataset(config, 'verbose', false)  % Override options
%
% Parameters:
%   config - Configuration object from create_sstdr_dataset_config()
%   varargin - Optional parameter overrides:
%             'verbose' - Display progress messages (default: true)
%             'save_progress' - Save progress every N networks (default: 10)

%% Parse inputs
if nargin < 1 || isempty(config)
    config = create_sstdr_dataset_config();
end

p = inputParser;
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'save_progress', 10, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});
opts = p.Results;

%% Add paths
addpath('functions/network_generation');
addpath('functions/sstdr_simulation');
addpath('config');

%% Disable Simulink warnings and GUI
warning('off', 'Simulink:Engine:MdlFileShadowedByFile');
warning('off', 'Simulink:Commands:LoadingOlderModel');
set_param(0, 'ShowLineWidths', 'off');
set_param(0, 'HideAutomaticNames', 'on');

%% Create output directory
if ~exist(config.dataset.output_dir, 'dir')
    mkdir(config.dataset.output_dir);
    if opts.verbose
        fprintf('Created output directory: %s\n', config.dataset.output_dir);
    end
end

%% Generate network configurations
if opts.verbose
    fprintf('\n--- Generating Network Configurations ---\n');
end

network_configs = cell(config.dataset.num_networks, 1);

for i = 1:config.dataset.num_networks
    % Generate random number of segments
    num_segments = randi(config.network.num_segments_range);
    
    % Generate network configuration
    network_configs{i} = generate_random_network( ...
        'num_segments', num_segments, ...
        'fault_probability', config.network.fault_probability, ...
        'fault_magnitude_range', config.network.fault_magnitude_range, ...
        'series_bias', config.network.series_bias, ...
        'dx', config.network.dx, ...
        'Z0', config.network.Z0, ...
        'velocity', config.network.velocity, ...
        'seed', i);
    
    % Add transmission line file to network config
    network_configs{i}.transmission_line_file = config.network.transmission_line_file;
    
    if opts.verbose && mod(i, 10) == 0
        fprintf('  Generated %d/%d network configurations\n', i, config.dataset.num_networks);
    end
end

if opts.verbose
    fprintf('Network configuration generation complete!\n');
end

%% Generate SSTDR data for each network
if opts.verbose
    fprintf('\n--- Generating SSTDR Data ---\n');
end

% Initialize results storage
dataset_results = struct();
dataset_results.config = config;
dataset_results.networks = cell(config.dataset.num_networks, 1);
dataset_results.metadata = struct();
dataset_results.metadata.created = datestr(now);
dataset_results.metadata.total_networks = config.dataset.num_networks;

% Performance tracking
start_time = tic;
processing_times = zeros(config.dataset.num_networks, 1);

for i = 1:config.dataset.num_networks
    network_start_time = tic;
    
    if opts.verbose
        fprintf('Processing network %d/%d...', i, config.dataset.num_networks);
    end
    
    net_config = network_configs{i};
    
    % Build temporary model
    model_name = sprintf('temp_network_%03d', i);
    
    try
        % Build network model (without opening GUI)
        build_network_model(net_config, ...
            'model_name', model_name, ...
            'save_model', false, ...
            'close_after', false, ...
            'connect_blocks', true);
        
        % Create SSTDR configuration from dataset config
        sstdr_cfg = create_sstdr_config_from_dataset_config(config);
        
        % Configure model sampling with exact max_step
        configure_model_sampling(model_name, 'max_step', config.simulation.max_step);
        
        % Set simulation to run without opening scope windows
        set_param(model_name, 'SimulationCommand', 'start');
        set_param(model_name, 'StopTime', num2str(config.simulation.duration));
        
        % Run simulation and analysis (skip plotting)
        run_sstdr_analysis('skip', model_name, 'save_results', false, 'plot_results', false);
        
        % Get results from base workspace
        sim_out = evalin('base', 'out');
        corr_results = evalin('base', 'correlation_results');
        
        % Store results
        network_result = struct();
        network_result.net_config = net_config;
        network_result.sstdr_cfg = sstdr_cfg;
        network_result.sim_out = sim_out;
        network_result.corr_results = corr_results;
        
        dataset_results.networks{i} = network_result;
        
        % Save individual file if requested
        if config.dataset.save_individual
            save_path = fullfile(config.dataset.output_dir, sprintf('network_%03d.mat', i));
            save(save_path, 'net_config', 'sstdr_cfg', 'sim_out', 'corr_results');
        end
        
        % Clean up model (close without saving)
        close_system(model_name, 0);
        
        processing_times(i) = toc(network_start_time);
        
        if opts.verbose
            fprintf(' %.2f s\n', processing_times(i));
        end
        
    catch ME
        if opts.verbose
            fprintf(' FAILED: %s\n', ME.message);
        end
        
        % Clean up model on error
        try
            close_system(model_name, 0);
        catch
            % Model may not exist
        end
        
        % Store error information
        network_result = struct();
        network_result.net_config = net_config;
        network_result.error = ME.message;
        dataset_results.networks{i} = network_result;
        
        processing_times(i) = toc(network_start_time);
    end
    
    % Save progress periodically
    if mod(i, opts.save_progress) == 0
        progress_file = fullfile(config.dataset.output_dir, sprintf('progress_%03d.mat', i));
        save(progress_file, 'dataset_results', 'processing_times');
        if opts.verbose
            fprintf('  Progress saved: %s\n', progress_file);
        end
    end
end

%% Calculate statistics and save summary
total_time = toc(start_time);
dataset_results.metadata.total_time = total_time;
dataset_results.metadata.avg_time_per_network = mean(processing_times);
dataset_results.metadata.processing_times = processing_times;

% Count successful networks
successful_networks = 0;
for i = 1:config.dataset.num_networks
    if isfield(dataset_results.networks{i}, 'corr_results')
        successful_networks = successful_networks + 1;
    end
end

dataset_results.metadata.successful_networks = successful_networks;
dataset_results.metadata.success_rate = successful_networks / config.dataset.num_networks;

%% Save final results
if config.dataset.save_summary
    summary_file = fullfile(config.dataset.output_dir, 'dataset_summary.mat');
    save(summary_file, 'dataset_results', 'config');
    if opts.verbose
        fprintf('\nDataset summary saved: %s\n', summary_file);
    end
end

%% Display final statistics
if opts.verbose
    fprintf('\n=== Dataset Generation Complete ===\n');
    fprintf('Configuration: %s\n', config.metadata.name);
    fprintf('Total Networks: %d\n', config.dataset.num_networks);
    fprintf('Successful: %d (%.1f%%)\n', successful_networks, dataset_results.metadata.success_rate * 100);
    fprintf('Total Time: %.1f minutes\n', total_time / 60);
    fprintf('Average Time per Network: %.2f seconds\n', dataset_results.metadata.avg_time_per_network);
    fprintf('Output Directory: %s\n', config.dataset.output_dir);
    
    % Display configuration summary
    fprintf('\nConfiguration Summary:\n');
    fprintf('  SSTDR: %.1f MHz chip rate, %.1f MHz sampling\n', config.sstdr.chip_rate/1e6, config.sstdr.fs/1e6);
    fprintf('  PN Sequence: %d bits (%d chips)\n', config.sstdr.pn_bits, 2^config.sstdr.pn_bits - 1);
    fprintf('  Simulation: %.1f μs duration\n', config.simulation.duration * 1e6);
    fprintf('  Networks: %d-%d segments, %.1f%% fault probability\n', ...
        config.network.num_segments_range(1), config.network.num_segments_range(2), ...
        config.network.fault_probability * 100);
    
    fprintf('\nDataset generation complete!\n');
end

end

function sstdr_cfg = create_sstdr_config_from_dataset_config(dataset_config)
%CREATE_SSTDR_CONFIG_FROM_DATASET_CONFIG Convert dataset config to SSTDR config

% Create SSTDR configuration using custom config
sstdr_cfg = sstdr_custom_config( ...
    'chip_rate', dataset_config.sstdr.chip_rate, ...
    'carrier_freq', dataset_config.sstdr.carrier_freq, ...
    'fs', dataset_config.sstdr.fs, ...
    'modulation', dataset_config.sstdr.modulation, ...
    'pn_bits', dataset_config.sstdr.pn_bits, ...
    'duration', dataset_config.simulation.duration, ...
    'apply', true);

% Override correlation settings from dataset config
if isfield(dataset_config, 'correlation')
    sstdr_cfg.correlation_config.method = dataset_config.correlation.method;
    sstdr_cfg.correlation_config.peak_threshold = dataset_config.correlation.peak_threshold;
    sstdr_cfg.correlation_config.normalize = dataset_config.correlation.normalize;
    sstdr_cfg.correlation_config.plot_results = dataset_config.correlation.plot_results;
    sstdr_cfg.correlation_config.positive_only = dataset_config.correlation.positive_only;
end

end 