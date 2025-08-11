function [results, sim_data] = run_simulation(varargin)
%RUN_SIMULATION Run SSTDR simulation on a Simscape model
%
% This function uses the standard sstdr_config system to run simulations.
% It accepts either a dataset config (from create_sstdr_dataset_config) or
% a direct sstdr_config and runs the simulation.
%
% Usage:
%   results = run_simulation()                      % Use default config and model
%   results = run_simulation(config)                % Use provided config, default model
%   results = run_simulation(config, 'my_model')    % Use provided config and model
%   results = run_simulation('my_model')            % Use default config, specified model

%% Parse inputs
p = inputParser;
addOptional(p, 'config', [], @(x) isempty(x) || isstruct(x));
addOptional(p, 'model_name', 'sstdr_basic', @ischar);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});
opts = p.Results;

% Handle case where first argument is a string (model name)
if ischar(opts.config)
    opts.model_name = opts.config;
    opts.config = [];
end

% Create default config if none provided
if isempty(opts.config)
    if opts.verbose
        fprintf('No config provided, creating default configuration...\n');
    end
    opts.config = create_simulation_config();
end

config = opts.config;
model_name = opts.model_name;

%% Initialize
results = struct();
results.success = false;
results.timestamp = datestr(now);

if opts.verbose
    fprintf('\n=== SSTDR Simulation: %s ===\n', model_name);
end

%% Apply configuration to workspace
% The config should already be in sstdr_config format
% (either from sstdr_config/sstdr_custom_config or create_sstdr_dataset_config)
if opts.verbose
    fprintf('Applying configuration to workspace...\n');
end

% Generate PN code and set up all workspace variables
pn_result = gen_pn_code( ...
    'modulation', config.pn_config.modulation, ...
    'carrier_freq', config.pn_config.carrier_freq, ...
    'chip_rate', config.pn_config.chip_rate, ...
    'fs', config.pn_config.fs, ...
    'pn_bits', config.pn_config.pn_bits, ...
    'polynomial', config.pn_config.polynomial, ...
    'magnitude', config.pn_config.magnitude, ...
    'export_to_base', true);

% Set up simulation parameters in base workspace
assignin('base', 'sim_stop_time', config.simulation_config.stop_time);
assignin('base', 'sim_solver', config.simulation_config.solver);
assignin('base', 'sim_max_step', config.simulation_config.max_step);
assignin('base', 'sim_fs', config.pn_config.fs);
assignin('base', 'sim_ts', 1/config.pn_config.fs);
assignin('base', 'sim_decimation', pn_result.settings.interpFactor);
assignin('base', 'sim_chip_rate', config.pn_config.chip_rate);

%% Run simulation
if opts.verbose
    fprintf('Running simulation...\n');
end

% Load model
if ~bdIsLoaded(model_name)
    load_system(model_name);
end

% Apply simulation parameters
set_param(model_name, 'StopTime', num2str(config.simulation_config.stop_time));
set_param(model_name, 'Solver', config.simulation_config.solver);
set_param(model_name, 'MaxStep', num2str(config.simulation_config.max_step));

% Run simulation
sim_start = tic;
sim_out = sim(model_name);
sim_time = toc(sim_start);

% Store results in base workspace
assignin('base', 'out', sim_out);

if opts.verbose
    fprintf('✓ Simulation completed in %.2f seconds\n', sim_time);
end

%% Extract signals
% TX signal is the PN code
tx_signal = struct('Time', pn_result.time_interp, 'Data', pn_result.pn_interp_modulated);

% RX signal from simulation - grab first signal and extract data
if isa(sim_out, 'Simulink.SimulationOutput')
    signal_names = sim_out.who;
    raw_signal = sim_out.get(signal_names{1});
else
    fields = fieldnames(sim_out);
    raw_signal = sim_out.(fields{1});
end

% Extract time and data from Simulink format
rx_signal = struct('Time', raw_signal.time, 'Data', raw_signal.signals.values);

if opts.verbose
    fprintf('✓ TX signal: %d samples\n', length(tx_signal.Data));
    fprintf('✓ RX signal: %d samples\n', length(rx_signal.Data));
end

%% Cross-correlation using the proper cross_correlate function
if opts.verbose
    fprintf('Running cross-correlation analysis...\n');
end

% Use the correlation_config from the passed config
if isfield(config, 'correlation_config')
    corr_config = config.correlation_config;
else
    % Default correlation configuration
    corr_config = struct('method', 'freq', 'plot_results', true, ...
                        'find_peaks', true, 'peak_threshold', 0.1, ...
                        'normalize', true, 'export_results', false, ...
                        'positive_only', false);
end

% Convert correlation config to name-value pairs
corr_params = {};
corr_fields = fieldnames(corr_config);
for i = 1:length(corr_fields)
    corr_params{end+1} = corr_fields{i};
    corr_params{end+1} = corr_config.(corr_fields{i});
end

% Run cross-correlation analysis
% Note: cross_correlate will get simulation output from base workspace 'out' variable
correlation_results = cross_correlate(corr_params{:});

if opts.verbose
    fprintf('✓ Cross-correlation analysis completed\n');
    
    % Display peak information based on method used
    if isfield(correlation_results, 'peak_times_time') && ~isempty(correlation_results.peak_times_time)
        fprintf('  Time domain peaks at: ');
        fprintf('%.3f ', correlation_results.peak_times_time * 1e6);
        fprintf('μs\n');
    end
    
    if isfield(correlation_results, 'peak_times_freq') && ~isempty(correlation_results.peak_times_freq)
        fprintf('  Frequency domain peaks at: ');
        fprintf('%.3f ', correlation_results.peak_times_freq * 1e6);
        fprintf('μs\n');
    end
end

%% Package results
results.success = true;
results.sim_time = sim_time;
results.tx_signal = tx_signal;
results.rx_signal = rx_signal;
results.correlation_results = correlation_results;

% Extract main correlation data based on method used
if isfield(correlation_results, 'correlation_time')
    results.correlation = correlation_results.correlation_time;
    results.time_axis = correlation_results.time_lags_time;
    if isfield(correlation_results, 'peaks_time')
        results.peaks = struct('values', correlation_results.peaks_time, ...
                              'times', correlation_results.peak_times_time);
    end
elseif isfield(correlation_results, 'correlation_freq')
    results.correlation = real(correlation_results.correlation_freq);
    results.time_axis = correlation_results.time_lags_freq;
    if isfield(correlation_results, 'peaks_freq')
        results.peaks = struct('values', correlation_results.peaks_freq, ...
                              'times', correlation_results.peak_times_freq);
    end
end

sim_data = sim_out;

end 