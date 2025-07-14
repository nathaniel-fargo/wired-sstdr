function config = create_sstdr_dataset_config(varargin)
%CREATE_SSTDR_DATASET_CONFIG Create configuration for SSTDR dataset generation
%
% This function extends sstdr_config with dataset-specific parameters.
% It uses the existing sstdr_config system as the base and adds network
% and dataset generation parameters on top.
%
% Usage:
%   config = create_sstdr_dataset_config()                    % Use defaults
%   config = create_sstdr_dataset_config('chip_rate', 10e6)   % Override SSTDR params
%   config = create_sstdr_dataset_config('num_networks', 50)  % Override dataset params
%
% SSTDR Parameters (passed to sstdr_config):
%   'chip_rate', 'carrier_freq', 'fs', 'pn_bits', 'modulation'
%
% Dataset Parameters:
%   'num_networks' - Number of networks to generate (default: 10)
%   'output_dir' - Output directory (default: auto-generated)
%   'dx' - Segment length in meters (default: 1000)
%   'velocity' - Propagation velocity in m/s (default: 2.75e8)
%   'fault_probability' - Fault probability per segment (default: 0.3)
%   'num_segments_range' - Range of segments per network (default: [5, 20])

%% Parse inputs
p = inputParser;

% SSTDR parameters (will be passed to sstdr_custom_config)
addParameter(p, 'chip_rate', 5e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'carrier_freq', 5e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 20e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'pn_bits', 11, @(x) isnumeric(x) && x >= 3 && x <= 20);
addParameter(p, 'modulation', 'sine', @(x) ismember(x, {'none', 'sine', 'cosine'}));
addParameter(p, 'duration', [], @(x) isempty(x) || (isnumeric(x) && x > 0));

% Dataset-specific parameters
addParameter(p, 'num_networks', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'output_dir', '', @ischar);

% Network parameters
addParameter(p, 'dx', 1000, @(x) isnumeric(x) && x > 0);
addParameter(p, 'velocity', 2.75e8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fault_probability', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'fault_magnitude_range', [0.1, 0.5], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'num_segments_range', [5, 20], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'transmission_line_file', 'line_test.mat', @ischar);

parse(p, varargin{:});
opts = p.Results;

%% Create base SSTDR configuration
fprintf('=== Creating SSTDR Dataset Configuration ===\n');

% Use sstdr_custom_config to create the base SSTDR config
sstdr_params = {
    'chip_rate', opts.chip_rate, ...
    'carrier_freq', opts.carrier_freq, ...
    'fs', opts.fs, ...
    'pn_bits', opts.pn_bits, ...
    'modulation', opts.modulation, ...
    'apply', false  % Don't apply to workspace yet
};

if ~isempty(opts.duration)
    sstdr_params = [sstdr_params, {'duration', opts.duration}];
end

% Create base SSTDR configuration
base_config = sstdr_custom_config(sstdr_params{:});

%% Add dataset-specific parameters
config = base_config;  % Start with base SSTDR config

% Add dataset metadata
config.dataset = struct();
config.dataset.num_networks = opts.num_networks;
config.dataset.save_individual = false;
config.dataset.save_summary = true;

% Auto-generate output directory if not provided
if isempty(opts.output_dir)
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    config.dataset.output_dir = fullfile('datasets', sprintf('sstdr_dataset_%s', timestamp));
else
    if ~startsWith(opts.output_dir, 'datasets/')
        config.dataset.output_dir = fullfile('datasets', opts.output_dir);
    else
        config.dataset.output_dir = opts.output_dir;
    end
end

% Add network parameters
config.network = struct();
config.network.dx = opts.dx;
config.network.velocity = opts.velocity;
config.network.Z0 = 50;  % Fixed impedance
config.network.fault_probability = opts.fault_probability;
config.network.fault_magnitude_range = opts.fault_magnitude_range;
config.network.series_bias = 0.5;  % 50/50 series/shunt
config.network.num_segments_range = opts.num_segments_range;
config.network.transmission_line_file = opts.transmission_line_file;

% Add correlation parameters (extend base config)
config.correlation_config.positive_only = true;  % Dataset generation uses positive only
config.correlation_config.plot_results = false;  % Don't plot during dataset generation

% Update metadata
config.name = sprintf('SSTDR Dataset Config (%.1f MHz)', opts.chip_rate/1e6);
config.type = 'dataset';  % Mark as dataset config

%% Display summary
fprintf('SSTDR Settings:\n');
fprintf('  Chip Rate: %.1f MHz\n', config.pn_config.chip_rate/1e6);
fprintf('  Carrier Freq: %.1f MHz\n', config.pn_config.carrier_freq/1e6);
fprintf('  Sampling Rate: %.1f MHz\n', config.pn_config.fs/1e6);
fprintf('  PN Bits: %d (%d chips)\n', config.pn_config.pn_bits, 2^config.pn_config.pn_bits-1);
fprintf('  Max Step: %.2f ns\n', config.simulation_config.max_step*1e9);

fprintf('\nNetwork Settings:\n');
fprintf('  Segment Length: %.1f km\n', config.network.dx/1000);
fprintf('  Velocity: %.2e m/s\n', config.network.velocity);
fprintf('  Segments per Network: %d-%d\n', config.network.num_segments_range(1), config.network.num_segments_range(2));
fprintf('  Fault Probability: %.1f%%\n', config.network.fault_probability*100);
fprintf('  Transmission Line: %s\n', config.network.transmission_line_file);

fprintf('\nDataset Settings:\n');
fprintf('  Number of Networks: %d\n', config.dataset.num_networks);
fprintf('  Output Directory: %s\n', config.dataset.output_dir);
fprintf('  Simulation Duration: %.1f μs\n', config.simulation_config.stop_time*1e6);

end 