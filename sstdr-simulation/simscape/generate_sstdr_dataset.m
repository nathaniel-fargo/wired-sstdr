% generate_sstdr_dataset - Generate SSTDR dataset from random networks
% 
% Usage:
%   generate_sstdr_dataset(num_networks, output_dir, sstdr_config_name)
%
% Inputs:
%   num_networks      - Number of networks to generate (default: 10)
%   output_dir        - Output directory name (default: 'dataset_<timestamp>')
%   sstdr_config_name - SSTDR configuration name (default: 'default')

function generate_sstdr_dataset(num_networks, output_dir, sstdr_config_name)

if nargin < 1 || isempty(num_networks)
    num_networks = 10;
end

if nargin < 2 || isempty(output_dir)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    output_dir = sprintf('dataset_%s', timestamp);
end

if nargin < 3 || isempty(sstdr_config_name)
    sstdr_config_name = 'default';
end

% Add paths
addpath('functions');
addpath('config');

fprintf('Generating dataset of %d networks in %s using %s config...\n', ...
    num_networks, output_dir, sstdr_config_name);

% Generate network configurations
network_configs = generate_network_dataset(num_networks);

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:num_networks
    fprintf('Processing network %d/%d...\n', i, num_networks);
    
    net_config = network_configs{i};
    
    % Build temporary model
    model_name = sprintf('temp_network_%03d', i);
    build_network_model(net_config, ...
        'model_name', model_name, ...
        'save_model', false, ...
        'close_after', false, ...
        'connect_blocks', true);
    
    % Load SSTDR configuration
    sstdr_cfg = sstdr_config(sstdr_config_name);
    
    % Configure model sampling
    configure_model_sampling(model_name);
    
    % Run simulation
    sim_out = sim(model_name);
    
    % Perform cross-correlation
    corr_results = cross_correlate(...
        'sim_output', sim_out, ...
        'method', 'freq', ...
        'plot_results', false);
    
    % Save data
    save_path = fullfile(output_dir, sprintf('network_%03d.mat', i));
    save(save_path, 'net_config', 'sstdr_cfg', 'sim_out', 'corr_results');
    
    fprintf('  Saved: %s\n', save_path);
    
    % Clean up model
    close_system(model_name, 0);
end

fprintf('Dataset generation complete! Output in: %s\n', output_dir);

end 