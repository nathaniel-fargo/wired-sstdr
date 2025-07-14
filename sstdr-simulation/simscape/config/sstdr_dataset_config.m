function [config] = sstdr_dataset_config(config_name)
%SSTDR_DATASET_CONFIG Configuration for SSTDR dataset generation
%
% Usage:
%   config = sstdr_dataset_config('high_freq_10mhz')  % Load predefined config
%   config = sstdr_dataset_config('custom')           % Load custom config
%
% This function provides easy-to-edit configurations for SSTDR dataset generation
% with all key parameters in one place.

if nargin < 1
    config_name = 'high_freq_10mhz';
end

fprintf('Loading SSTDR dataset configuration: %s\n', config_name);

%% Configuration Selection
switch lower(config_name)
    case 'high_freq_10mhz'
        % 10 MHz matched chip/carrier rate configuration
        config = create_high_freq_10mhz_config();
        
    case 'custom'
        % Custom configuration - edit this section for your needs
        config = create_custom_config();
        
    case 'legacy_compatible'
        % Legacy compatible configuration for comparison
        config = create_legacy_config();
        
    otherwise
        error('Unknown configuration: %s. Available: high_freq_10mhz, custom, legacy_compatible', config_name);
end

%% Display Configuration Summary
display_config_summary(config);

end

function config = create_high_freq_10mhz_config()
%CREATE_HIGH_FREQ_10MHZ_CONFIG 10 MHz matched chip/carrier rate configuration

config = struct();

%% SSTDR Signal Parameters
config.sstdr = struct();
config.sstdr.chip_rate = 10e6;           % 10 MHz chip rate
config.sstdr.carrier_freq = 10e6;        % 10 MHz carrier (matched to chip rate)
config.sstdr.fs = 40e6;                  % 40 MHz sampling (4x chip rate)
config.sstdr.modulation = 'sine';        % Sine wave modulation
config.sstdr.pn_bits = 10;               % 2^10 - 1 = 1023 chips
config.sstdr.polynomial = [10 3 0];      % 10-bit primitive polynomial
config.sstdr.magnitude = 1;              % Signal magnitude

%% Simulation Parameters
config.simulation = struct();
config.simulation.duration = 200e-6;     % 200 microseconds simulation time
config.simulation.solver = 'ode45';      % Solver for smooth signals
config.simulation.max_step = 1/(10*config.sstdr.fs);  % 10 steps per sample

%% Cross-Correlation Parameters
config.correlation = struct();
config.correlation.method = 'freq';      % Use frequency domain for speed
config.correlation.peak_threshold = 0.1; % Peak detection threshold
config.correlation.normalize = true;     % Normalize correlation output
config.correlation.positive_only = true; % Return only 0 to +t (NEW FEATURE)
config.correlation.plot_results = false; % Don't plot during dataset generation

%% Network Generation Parameters
config.network = struct();
config.network.num_segments_range = [5, 15];     % Random segments between 5-15
config.network.fault_probability = 0.2;          % 20% chance of fault per segment
config.network.fault_magnitude_range = [0.1, 1.0]; % Fault magnitude range
config.network.series_bias = 0.5;                % 50% series, 50% shunt faults
config.network.dx = 1.0;                         % 1 meter per segment
config.network.Z0 = 50;                          % 50 ohm characteristic impedance
config.network.velocity = 2e8;                   % 2e8 m/s propagation velocity
config.network.transmission_line_file = 'line_test.mat'; % Transmission line model file

%% Dataset Generation Parameters
config.dataset = struct();
config.dataset.num_networks = 10;                % Number of networks to generate
config.dataset.output_dir = 'dataset_10mhz';     % Output directory
config.dataset.save_individual = true;           % Save individual network files
config.dataset.save_summary = true;              % Save dataset summary

%% Configuration Metadata
config.metadata = struct();
config.metadata.name = 'High Frequency 10 MHz SSTDR';
config.metadata.description = '10 MHz matched chip/carrier rate with 40 MHz sampling';
config.metadata.created = datestr(now);
config.metadata.version = '1.0';

end

function config = create_custom_config()
%CREATE_CUSTOM_CONFIG Custom configuration - EDIT THIS FOR YOUR NEEDS

config = struct();

%% SSTDR Signal Parameters - EDIT THESE VALUES
config.sstdr = struct();
config.sstdr.chip_rate = 5e6;            % EDIT: Chip rate in Hz
config.sstdr.carrier_freq = 5e6;         % EDIT: Carrier frequency in Hz
config.sstdr.fs = 20e6;                  % EDIT: Sampling frequency in Hz
config.sstdr.modulation = 'sine';        % EDIT: 'sine', 'cosine', or 'none'
config.sstdr.pn_bits = 11;               % EDIT: PN sequence bits (2^11-1 = 2047 chips)
config.sstdr.polynomial = [11 2 0];      % EDIT: Primitive polynomial for PN sequence
config.sstdr.magnitude = 1;              % EDIT: Signal magnitude

%% Simulation Parameters - EDIT THESE VALUES
config.simulation = struct();
config.simulation.duration = 500e-6;     % EDIT: Simulation duration in seconds
config.simulation.solver = 'ode45';      % EDIT: Solver type
config.simulation.max_step = 1/(10*config.sstdr.fs);

%% Cross-Correlation Parameters - EDIT THESE VALUES
config.correlation = struct();
config.correlation.method = 'freq';      % EDIT: 'time', 'freq', or 'both'
config.correlation.peak_threshold = 0.15; % EDIT: Peak detection threshold
config.correlation.normalize = true;     % EDIT: Normalize correlation output
config.correlation.positive_only = true; % EDIT: Return only 0 to +t
config.correlation.plot_results = false; % EDIT: Plot results during generation

%% Network Generation Parameters - EDIT THESE VALUES
config.network = struct();
config.network.num_segments_range = [8, 20];     % EDIT: Random segments range
config.network.fault_probability = 0.25;         % EDIT: Fault probability per segment
config.network.fault_magnitude_range = [0.2, 0.8]; % EDIT: Fault magnitude range
config.network.series_bias = 0.6;                % EDIT: Series vs shunt bias
config.network.dx = 0.5;                         % EDIT: Meters per segment
config.network.Z0 = 75;                          % EDIT: Characteristic impedance
config.network.velocity = 2.2e8;                 % EDIT: Propagation velocity
config.network.transmission_line_file = 'line_test.mat'; % EDIT: T-line model file

%% Dataset Generation Parameters - EDIT THESE VALUES
config.dataset = struct();
config.dataset.num_networks = 50;                % EDIT: Number of networks
config.dataset.output_dir = 'dataset_custom';    % EDIT: Output directory
config.dataset.save_individual = true;           % EDIT: Save individual files
config.dataset.save_summary = true;              % EDIT: Save summary

%% Configuration Metadata
config.metadata = struct();
config.metadata.name = 'Custom SSTDR Configuration';
config.metadata.description = 'User-customized SSTDR dataset configuration';
config.metadata.created = datestr(now);
config.metadata.version = '1.0';

end

function config = create_legacy_config()
%CREATE_LEGACY_CONFIG Legacy compatible configuration for comparison

config = struct();

%% SSTDR Signal Parameters (legacy defaults)
config.sstdr = struct();
config.sstdr.chip_rate = 100e3;          % 100 kHz chip rate
config.sstdr.carrier_freq = 100e3;       % 100 kHz carrier
config.sstdr.fs = 400e3;                 % 400 kHz sampling
config.sstdr.modulation = 'sine';        % Sine wave modulation
config.sstdr.pn_bits = 10;               % 1023 chips
config.sstdr.polynomial = [10 3 0];      % 10-bit primitive polynomial
config.sstdr.magnitude = 1;              % Signal magnitude

%% Simulation Parameters
config.simulation = struct();
config.simulation.duration = 20e-3;      % 20 milliseconds (longer for low freq)
config.simulation.solver = 'ode45';      % Solver
config.simulation.max_step = 1/(10*config.sstdr.fs);

%% Cross-Correlation Parameters
config.correlation = struct();
config.correlation.method = 'freq';      % Frequency domain
config.correlation.peak_threshold = 0.1; % Peak detection threshold
config.correlation.normalize = true;     % Normalize correlation output
config.correlation.positive_only = false; % Keep legacy -t to +t behavior
config.correlation.plot_results = false; % Don't plot during generation

%% Network Generation Parameters
config.network = struct();
config.network.num_segments_range = [5, 15];     % Random segments
config.network.fault_probability = 0.15;         % 15% fault probability
config.network.fault_magnitude_range = [0.1, 1.0]; % Fault magnitude range
config.network.series_bias = 0.5;                % Equal series/shunt
config.network.dx = 1.0;                         % 1 meter per segment
config.network.Z0 = 50;                          % 50 ohm impedance
config.network.velocity = 2e8;                   % 2e8 m/s velocity
config.network.transmission_line_file = 'line_test.mat'; % T-line model

%% Dataset Generation Parameters
config.dataset = struct();
config.dataset.num_networks = 10;                % Number of networks
config.dataset.output_dir = 'dataset_legacy';    % Output directory
config.dataset.save_individual = true;           % Save individual files
config.dataset.save_summary = true;              % Save summary

%% Configuration Metadata
config.metadata = struct();
config.metadata.name = 'Legacy Compatible SSTDR';
config.metadata.description = 'Legacy compatible configuration for comparison';
config.metadata.created = datestr(now);
config.metadata.version = '1.0';

end

function display_config_summary(config)
%DISPLAY_CONFIG_SUMMARY Display configuration summary

fprintf('\n=== SSTDR Dataset Configuration Summary ===\n');
fprintf('Name: %s\n', config.metadata.name);
fprintf('Description: %s\n', config.metadata.description);

fprintf('\nSSTDR Signal Parameters:\n');
fprintf('  Chip Rate: %.1f MHz (%.3f μs per chip)\n', config.sstdr.chip_rate/1e6, 1e6/config.sstdr.chip_rate);
fprintf('  Carrier Frequency: %.1f MHz\n', config.sstdr.carrier_freq/1e6);
fprintf('  Sampling Frequency: %.1f MHz\n', config.sstdr.fs/1e6);
fprintf('  Modulation: %s\n', config.sstdr.modulation);
fprintf('  PN Bits: %d (sequence length: %d chips)\n', config.sstdr.pn_bits, 2^config.sstdr.pn_bits - 1);
fprintf('  PN Duration: %.1f μs\n', (2^config.sstdr.pn_bits - 1) / config.sstdr.chip_rate * 1e6);

fprintf('\nSimulation Parameters:\n');
fprintf('  Duration: %.1f μs\n', config.simulation.duration * 1e6);
fprintf('  Solver: %s\n', config.simulation.solver);
fprintf('  Max Step: %.1f ns\n', config.simulation.max_step * 1e9);

fprintf('\nCorrelation Parameters:\n');
fprintf('  Method: %s domain\n', config.correlation.method);
fprintf('  Peak Threshold: %.2f\n', config.correlation.peak_threshold);
fprintf('  Normalize: %s\n', bool2str(config.correlation.normalize));
fprintf('  Positive Only: %s\n', bool2str(config.correlation.positive_only));

fprintf('\nNetwork Parameters:\n');
fprintf('  Segments: %d-%d (random)\n', config.network.num_segments_range(1), config.network.num_segments_range(2));
fprintf('  Fault Probability: %.1f%%\n', config.network.fault_probability * 100);
fprintf('  Fault Magnitude: %.1f-%.1f\n', config.network.fault_magnitude_range(1), config.network.fault_magnitude_range(2));
fprintf('  Series Bias: %.1f%% series, %.1f%% shunt\n', config.network.series_bias * 100, (1-config.network.series_bias) * 100);
fprintf('  Segment Length: %.1f m\n', config.network.dx);
fprintf('  Characteristic Impedance: %.1f Ω\n', config.network.Z0);
fprintf('  Propagation Velocity: %.1e m/s\n', config.network.velocity);
fprintf('  Transmission Line File: %s\n', config.network.transmission_line_file);

fprintf('\nDataset Parameters:\n');
fprintf('  Number of Networks: %d\n', config.dataset.num_networks);
fprintf('  Output Directory: %s\n', config.dataset.output_dir);
fprintf('  Save Individual Files: %s\n', bool2str(config.dataset.save_individual));
fprintf('  Save Summary: %s\n', bool2str(config.dataset.save_summary));

fprintf('\nFrequency Relationships:\n');
fprintf('  Sampling Ratio: %.1fx (fs/chip_rate)\n', config.sstdr.fs/config.sstdr.chip_rate);
fprintf('  Carrier/Chip Ratio: %.1fx\n', config.sstdr.carrier_freq/config.sstdr.chip_rate);

end

function str = bool2str(val)
%BOOL2STR Convert boolean to string
if val
    str = 'Yes';
else
    str = 'No';
end
end 