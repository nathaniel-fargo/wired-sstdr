function config = create_sstdr_dataset_config(varargin)
%CREATE_SSTDR_DATASET_CONFIG Create configuration object for SSTDR dataset generation
%
% Usage:
%   config = create_sstdr_dataset_config()                    % Use defaults
%   config = create_sstdr_dataset_config('chip_rate', 10e6)   % Override parameters
%
% Default Parameters:
%   chip_rate: 5 MHz
%   carrier_freq: 5 MHz  
%   fs: 20 MHz
%   pn_bits: 11 (2047 chips)
%   max_step: 1/fs
%   dx: 1 km
%   velocity: 2.75e8 m/s
%
% Optional Parameters:
%   'chip_rate' - Chip rate in Hz (default: 5e6)
%   'carrier_freq' - Carrier frequency in Hz (default: 5e6)
%   'fs' - Sampling frequency in Hz (default: 20e6)
%   'pn_bits' - PN sequence bits (default: 11)
%   'duration' - Simulation duration in seconds (default: auto from PN)
%   'dx' - Segment length in meters (default: 1000)
%   'velocity' - Propagation velocity in m/s (default: 2.75e8)
%   'num_networks' - Number of networks to generate (default: 10)
%   'output_dir' - Output directory (default: 'sstdr_dataset')
%   'fault_probability' - Fault probability per segment (default: 0.3)
%   'fault_magnitude_range' - Fault magnitude range (default: [0.1, 0.5])
%   'num_segments_range' - Range of segments per network (default: [5, 20])
%   'transmission_line_file' - Transmission line file (default: 'line_test.mat')

%% Parse inputs
p = inputParser;

% SSTDR parameters
addParameter(p, 'chip_rate', 5e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'carrier_freq', 5e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 20e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'pn_bits', 11, @(x) isnumeric(x) && x >= 3 && x <= 20);
addParameter(p, 'duration', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'modulation', 'bpsk', @ischar);

% Network parameters
addParameter(p, 'dx', 1000, @(x) isnumeric(x) && x > 0);
addParameter(p, 'velocity', 2.75e8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'Z0', 50, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fault_probability', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'fault_magnitude_range', [0.1, 0.5], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'series_bias', 0.5, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'num_segments_range', [5, 20], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'transmission_line_file', 'line_test.mat', @ischar);

% Dataset parameters
addParameter(p, 'num_networks', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'output_dir', 'sstdr_dataset', @ischar);
addParameter(p, 'save_individual', false, @islogical);
addParameter(p, 'save_summary', true, @islogical);

% Correlation parameters
addParameter(p, 'correlation_method', 'freq', @ischar);
addParameter(p, 'positive_only', true, @islogical);
addParameter(p, 'normalize', true, @islogical);
addParameter(p, 'peak_threshold', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);

parse(p, varargin{:});
opts = p.Results;

%% Calculate derived parameters
% Auto-calculate duration if not provided
if isempty(opts.duration)
    pn_length = 2^opts.pn_bits - 1;
    opts.duration = pn_length / opts.chip_rate * 3; % 3x PN sequence length
end

% Calculate max_step as exactly 1/fs
max_step = 1 / opts.fs;

%% Build configuration structure
config = struct();

% Metadata
config.metadata = struct();
config.metadata.name = 'SSTDR Dataset Configuration';
config.metadata.created = datestr(now);
config.metadata.version = '1.0';

% SSTDR configuration
config.sstdr = struct();
config.sstdr.chip_rate = opts.chip_rate;
config.sstdr.carrier_freq = opts.carrier_freq;
config.sstdr.fs = opts.fs;
config.sstdr.pn_bits = opts.pn_bits;
config.sstdr.modulation = opts.modulation;

% Simulation configuration
config.simulation = struct();
config.simulation.duration = opts.duration;
config.simulation.max_step = max_step;

% Network configuration
config.network = struct();
config.network.dx = opts.dx;
config.network.velocity = opts.velocity;
config.network.Z0 = opts.Z0;
config.network.fault_probability = opts.fault_probability;
config.network.fault_magnitude_range = opts.fault_magnitude_range;
config.network.series_bias = opts.series_bias;
config.network.num_segments_range = opts.num_segments_range;
config.network.transmission_line_file = opts.transmission_line_file;

% Dataset configuration
config.dataset = struct();
config.dataset.num_networks = opts.num_networks;
config.dataset.output_dir = opts.output_dir;
config.dataset.save_individual = opts.save_individual;
config.dataset.save_summary = opts.save_summary;

% Correlation configuration
config.correlation = struct();
config.correlation.method = opts.correlation_method;
config.correlation.positive_only = opts.positive_only;
config.correlation.normalize = opts.normalize;
config.correlation.peak_threshold = opts.peak_threshold;
config.correlation.plot_results = false; % Never plot during generation

%% Display configuration summary
fprintf('=== SSTDR Dataset Configuration ===\n');
fprintf('SSTDR Settings:\n');
fprintf('  Chip Rate: %.1f MHz\n', config.sstdr.chip_rate/1e6);
fprintf('  Carrier Freq: %.1f MHz\n', config.sstdr.carrier_freq/1e6);
fprintf('  Sampling Rate: %.1f MHz\n', config.sstdr.fs/1e6);
fprintf('  PN Bits: %d (%d chips)\n', config.sstdr.pn_bits, 2^config.sstdr.pn_bits - 1);
fprintf('  Max Step: %.2f ns\n', config.simulation.max_step * 1e9);

fprintf('\nNetwork Settings:\n');
fprintf('  Segment Length: %.1f km\n', config.network.dx/1000);
fprintf('  Velocity: %.2e m/s\n', config.network.velocity);
fprintf('  Segments per Network: %d-%d\n', config.network.num_segments_range(1), config.network.num_segments_range(2));
fprintf('  Fault Probability: %.1f%%\n', config.network.fault_probability * 100);

fprintf('\nDataset Settings:\n');
fprintf('  Number of Networks: %d\n', config.dataset.num_networks);
fprintf('  Output Directory: %s\n', config.dataset.output_dir);
fprintf('  Simulation Duration: %.1f μs\n', config.simulation.duration * 1e6);

end 