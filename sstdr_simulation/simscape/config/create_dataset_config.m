function dataset_config = create_dataset_config(varargin)
%CREATE_DATASET_CONFIG Create configuration for SSTDR dataset generation
%
% This function creates a configuration structure for generating datasets
% of random SSTDR networks with specified parameters for network size,
% fault characteristics, and simulation settings.
%
% Usage:
%   dataset_config = create_dataset_config()  % Use defaults
%   dataset_config = create_dataset_config('param', value, ...)
%
% Parameters:
%   name              - Dataset name (default: 'sstdr_dataset')
%   num_networks      - Number of networks to generate (default: 100)
%   
%   Network Generation Parameters:
%   segments_range    - [min, max] segments per network (default: [3, 10])
%   length_range      - [min, max] total length in meters (default: [10, 100])
%   fault_prob        - Probability of fault per segment (default: 0.3)
%   fault_magnitude_range - [min, max] fault magnitude (default: [0.1, 0.8])
%   fault_type_prob   - [open_prob, short_prob] (default: [0.6, 0.4])
%   termination_range - [min, max] termination impedance (default: [25, 100])
%   
%   Simulation Parameters:
%   network_config    - Base network config (optional)
%   simulation_config - Base simulation config (optional)
%   
%   Output Options:
%   data_folder       - Custom data folder name (default: auto-generated)
%   verbose           - Display progress (default: true)
%   save_models       - Save Simulink models (default: false)
%   parallel          - Use parallel processing (default: false)
%
% Output:
%   dataset_config - Configuration structure for dataset generation

%% Parse input arguments
p = inputParser;
addParameter(p, 'name', 'sstdr_dataset', @ischar);
addParameter(p, 'num_networks', 100, @(x) isnumeric(x) && x > 0);

% Network generation parameters
addParameter(p, 'segments_range', [3, 10], @(x) isnumeric(x) && length(x) == 2 && x(1) <= x(2));
addParameter(p, 'length_range', [10, 100], @(x) isnumeric(x) && length(x) == 2 && x(1) <= x(2));
addParameter(p, 'fault_prob', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'fault_magnitude_range', [0.1, 0.8], @(x) isnumeric(x) && length(x) == 2 && x(1) <= x(2));
addParameter(p, 'fault_type_prob', [0.6, 0.4], @(x) isnumeric(x) && length(x) == 2 && sum(x) == 1);
addParameter(p, 'termination_range', [25, 100], @(x) isnumeric(x) && length(x) == 2 && x(1) <= x(2));

% Simulation parameters
addParameter(p, 'network_config', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'simulation_config', [], @(x) isempty(x) || isstruct(x));

% Output options
addParameter(p, 'data_folder', '', @ischar);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'save_models', false, @islogical);
addParameter(p, 'parallel', false, @islogical);

parse(p, varargin{:});
opts = p.Results;

%% Create dataset configuration structure
dataset_config = struct();

% Basic info
dataset_config.name = opts.name;
dataset_config.num_networks = opts.num_networks;
dataset_config.created = datestr(now);

% Network generation parameters
dataset_config.network_generation = struct();
dataset_config.network_generation.segments_range = opts.segments_range;
dataset_config.network_generation.length_range = opts.length_range;
dataset_config.network_generation.fault_prob = opts.fault_prob;
dataset_config.network_generation.fault_magnitude_range = opts.fault_magnitude_range;
dataset_config.network_generation.fault_type_prob = opts.fault_type_prob;
dataset_config.network_generation.termination_range = opts.termination_range;

% Base configurations
if isempty(opts.network_config)
    % Create default network config template
    dataset_config.base_network_config = struct();
    dataset_config.base_network_config.type = 'random_generated';
    dataset_config.base_network_config.physical = struct();
    dataset_config.base_network_config.physical.Z0 = 50;  % Default characteristic impedance
    dataset_config.base_network_config.physical.velocity = 2e8;  % Default velocity
else
    dataset_config.base_network_config = opts.network_config;
end

if isempty(opts.simulation_config)
    % Create default simulation config with no plots
    dataset_config.base_simulation_config = create_simulation_config(...
        'name', 'dataset_generation', ...
        'plot_results', false, ...  % Disable plots for batch generation
        'apply', false);
else
    dataset_config.base_simulation_config = opts.simulation_config;
    % Force disable plots for batch generation
    dataset_config.base_simulation_config.correlation_config.plot_results = false;
end

% Output options
dataset_config.output = struct();
dataset_config.output.data_folder = opts.data_folder;
dataset_config.output.verbose = opts.verbose;
dataset_config.output.save_models = opts.save_models;
dataset_config.output.parallel = opts.parallel;

%% Display configuration summary
if opts.verbose
    fprintf('\n=== Dataset Generation Configuration ===\n');
    fprintf('Name: %s\n', dataset_config.name);
    fprintf('Networks to generate: %d\n', dataset_config.num_networks);
    fprintf('\nNetwork Parameters:\n');
    fprintf('  Segments: %d - %d\n', opts.segments_range(1), opts.segments_range(2));
    fprintf('  Length: %.1f - %.1f m\n', opts.length_range(1), opts.length_range(2));
    fprintf('  Fault probability: %.1f%% per segment\n', opts.fault_prob * 100);
    fprintf('  Fault magnitude: %.1f - %.1f\n', opts.fault_magnitude_range(1), opts.fault_magnitude_range(2));
    fprintf('  Fault types: %.1f%% open, %.1f%% short\n', opts.fault_type_prob(1)*100, opts.fault_type_prob(2)*100);
    fprintf('  Termination: %.1f - %.1f Ω\n', opts.termination_range(1), opts.termination_range(2));
    fprintf('\nOutput Options:\n');
    fprintf('  Verbose: %s\n', string(opts.verbose));
    fprintf('  Save models: %s\n', string(opts.save_models));
    fprintf('  Parallel processing: %s\n', string(opts.parallel));
    fprintf('  Plots disabled: true (for batch generation)\n');
end

end 