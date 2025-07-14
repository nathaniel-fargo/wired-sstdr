function [results, sim_data] = run_sstdr_simulation(config, model_name, varargin)
%RUN_SSTDR_SIMULATION Run SSTDR simulation on a Simscape model
%
% This function uses the standard sstdr_config system to run simulations.
% It accepts either a dataset config (from create_sstdr_dataset_config) or
% a direct sstdr_config and runs the simulation.
%
% Usage:
%   results = run_sstdr_simulation(config, 'my_model')
%   results = run_sstdr_simulation(config, 'my_model', 'plot_results', true)

%% Parse inputs
p = inputParser;
addRequired(p, 'config', @isstruct);
addRequired(p, 'model_name', @ischar);
addParameter(p, 'plot_results', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, config, model_name, varargin{:});
opts = p.Results;

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

%% Cross-correlation
if opts.verbose
    fprintf('Running cross-correlation...\n');
end

[correlation, lags] = xcorr(rx_signal.Data, tx_signal.Data);
time_axis = lags / config.pn_config.fs;

% Keep only positive delays (if configured)
if isfield(config, 'correlation_config') && isfield(config.correlation_config, 'positive_only') && config.correlation_config.positive_only
    positive_idx = time_axis >= 0;
    correlation = correlation(positive_idx);
    time_axis = time_axis(positive_idx);
end

% Normalize
correlation = correlation / max(abs(correlation));

% Find peaks
peak_threshold = 0.1;
if isfield(config, 'correlation_config') && isfield(config.correlation_config, 'peak_threshold')
    peak_threshold = config.correlation_config.peak_threshold;
end

[peak_vals, peak_locs] = findpeaks(abs(correlation), 'MinPeakHeight', peak_threshold);

if opts.verbose
    fprintf('✓ Found %d peaks\n', length(peak_vals));
end

%% Package results
results.success = true;
results.sim_time = sim_time;
results.tx_signal = tx_signal;
results.rx_signal = rx_signal;
results.correlation = correlation;
results.time_axis = time_axis;
results.peaks = struct('values', peak_vals, 'times', time_axis(peak_locs));

%% Plot if requested
if opts.plot_results
    figure('Name', 'SSTDR Results');
    
    subplot(2,1,1);
    plot(tx_signal.Time*1e6, tx_signal.Data, 'b-', rx_signal.Time*1e6, rx_signal.Data, 'r-');
    xlabel('Time (μs)'); ylabel('Amplitude');
    legend('TX', 'RX'); title('Signals');
    
    subplot(2,1,2);
    plot(time_axis*1e6, correlation, 'k-');
    hold on;
    if ~isempty(peak_vals)
        plot(time_axis(peak_locs)*1e6, peak_vals, 'ro', 'MarkerFaceColor', 'r');
    end
    xlabel('Time Delay (μs)'); ylabel('Correlation');
    title('Cross-Correlation'); grid on;
end

sim_data = sim_out;

end 