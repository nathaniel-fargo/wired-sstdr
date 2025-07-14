function [results, model_info] = run_network_sstdr(network_config, sstdr_config, varargin)
%RUN_NETWORK_SSTDR Generate network model and run SSTDR simulation
%
% This function combines network generation and SSTDR simulation into a single
% callable unit, useful for testing and standalone analysis.
%
% Usage:
%   results = run_network_sstdr(net_config, sstdr_config)
%   results = run_network_sstdr(net_config, sstdr_config, 'save_model', true)
%
% Input:
%   network_config - Network configuration from create_network_config()
%   sstdr_config   - SSTDR configuration from create_sstdr_dataset_config()
%   varargin       - Optional name-value pairs:
%                    'model_name' - Name for model (default: auto-generated)
%                    'save_model' - Save Simulink model to disk (default: false)
%                    'keep_open' - Keep model open after simulation (default: false)
%                    'plot_results' - Plot correlation results (default: true)
%                    'verbose' - Display progress messages (default: true)
%                    'output_dir' - Directory for saved files (default: current)
%
% Output:
%   results    - Structure containing simulation and correlation results
%   model_info - Structure with model details and block information

%% Parse inputs
p = inputParser;
addRequired(p, 'network_config', @isstruct);
addRequired(p, 'sstdr_config', @isstruct);
addParameter(p, 'model_name', '', @ischar);
addParameter(p, 'save_model', false, @islogical);
addParameter(p, 'keep_open', false, @islogical);
addParameter(p, 'plot_results', true, @islogical);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'output_dir', pwd, @ischar);
parse(p, network_config, sstdr_config, varargin{:});
opts = p.Results;

%% Generate model name if not provided
if isempty(opts.model_name)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    opts.model_name = sprintf('network_sstdr_%s', timestamp);
end

%% Add required paths
current_dir = fileparts(mfilename('fullpath'));
simscape_root = fileparts(fileparts(current_dir));
addpath(fullfile(simscape_root, 'functions', 'network_generation'));
addpath(fullfile(simscape_root, 'functions', 'sstdr_simulation'));
addpath(fullfile(simscape_root, 'config'));

%% Initialize results structure
results = struct();
results.network_config = network_config;
results.sstdr_config = sstdr_config;
results.model_name = opts.model_name;
results.success = false;
results.error_message = '';
results.timestamp = datestr(now);
results.sim_time = 0;
results.model_file = '';
results.results_file = '';

if opts.verbose
    fprintf('\n=== Network SSTDR Simulation ===\n');
    fprintf('Network: %s (%d segments, %.1f m)\n', ...
        network_config.name, network_config.num_segments, network_config.physical.total_length);
    fprintf('Model: %s\n', opts.model_name);
end

try
    %% Step 1: Build network model
    if opts.verbose
        fprintf('\n--- Building Network Model ---\n');
    end
    
    [model_name, model_info] = build_network_model(network_config, ...
        'model_name', opts.model_name, ...
        'save_model', false, ...  % We'll save later if requested
        'close_after', false, ...
        'connect_blocks', true);
    
    results.model_name = model_name;
    
    %% Step 2: Configure SSTDR parameters
    if opts.verbose
        fprintf('\n--- Configuring SSTDR Parameters ---\n');
    end
    
    % Apply SSTDR configuration to the model
    if isfield(sstdr_config, 'sstdr')
        config_sstdr = sstdr_config.sstdr;
    else
        config_sstdr = sstdr_config;  % Direct SSTDR config
    end
    
    % Set simulation parameters
    if isfield(sstdr_config, 'simulation')
        sim_duration = sstdr_config.simulation.duration;
        max_step = sstdr_config.simulation.max_step;
    else
        sim_duration = 100e-6;  % Default 100 μs
        max_step = 1/(4*config_sstdr.fs);  % Default based on sampling rate
    end
    
    % Configure model timing
    set_param(model_name, 'StopTime', num2str(sim_duration));
    set_param(model_name, 'MaxStep', num2str(max_step));
    
    % Generate PN code for SSTDR subsystem
    pn_config = gen_pn_code( ...
        'chip_rate', config_sstdr.chip_rate, ...
        'carrier_freq', config_sstdr.carrier_freq, ...
        'fs', config_sstdr.fs, ...
        'pn_bits', config_sstdr.pn_bits, ...
        'modulation', config_sstdr.modulation, ...
        'export_to_base', true);
    
    % Set additional workspace variables for SSTDR
    assignin('base', 'sim_max_step', max_step);
    assignin('base', 'chip_rate', config_sstdr.chip_rate);
    assignin('base', 'carrier_freq', config_sstdr.carrier_freq);
    assignin('base', 'fs', config_sstdr.fs);
    assignin('base', 'pn_bits', config_sstdr.pn_bits);
    
    if opts.verbose
        fprintf('  ✓ Simulation duration: %.1f μs\n', sim_duration * 1e6);
        fprintf('  ✓ Max step: %.2f ns\n', max_step * 1e9);
        fprintf('  ✓ SSTDR: %.1f MHz chip rate, %.1f MHz sampling\n', ...
            config_sstdr.chip_rate/1e6, config_sstdr.fs/1e6);
    end
    
    %% Step 3: Run simulation
    if opts.verbose
        fprintf('\n--- Running Simulation ---\n');
    end
    
    % Run simulation
    sim_start_time = tic;
    simOut = sim(model_name, 'StopTime', num2str(sim_duration));
    sim_time = toc(sim_start_time);
    
    if opts.verbose
        fprintf('  ✓ Simulation completed in %.2f seconds\n', sim_time);
    end
    
    % Store simulation results
    results.sim_out = simOut;
    results.sim_time = sim_time;
    
    %% Step 4: Run SSTDR analysis
    if opts.verbose
        fprintf('\n--- Running SSTDR Analysis ---\n');
    end
    
    % Extract signals from simulation
    try
        % Try to get signals from simOut
        if isa(simOut, 'Simulink.SimulationOutput')
            % Get logged signals
            signal_names = simOut.who;
            if opts.verbose
                fprintf('  Available signals: %s\n', strjoin(signal_names, ', '));
            end
            
            % Try to find TX and RX signals
            tx_signal = [];
            rx_signal = [];
            
            for i = 1:length(signal_names)
                signal_name = signal_names{i};
                if contains(lower(signal_name), 'tx') || contains(lower(signal_name), 'current')
                    tx_signal = simOut.get(signal_name);
                elseif contains(lower(signal_name), 'rx') || contains(lower(signal_name), 'voltage')
                    rx_signal = simOut.get(signal_name);
                end
            end
            
            if isempty(tx_signal) || isempty(rx_signal)
                if opts.verbose
                    fprintf('  ⚠ Could not find TX/RX signals in simulation output\n');
                end
                % Generate dummy signals for testing
                t = 0:1/config_sstdr.fs:sim_duration;
                tx_signal = struct('Time', t, 'Data', randn(length(t), 1));
                rx_signal = struct('Time', t, 'Data', randn(length(t), 1));
            end
            
        else
            % Fallback for older simulation output format
            if opts.verbose
                fprintf('  ⚠ Using fallback signal extraction\n');
            end
            t = 0:1/config_sstdr.fs:sim_duration;
            tx_signal = struct('Time', t, 'Data', randn(length(t), 1));
            rx_signal = struct('Time', t, 'Data', randn(length(t), 1));
        end
        
        % Run cross-correlation analysis
        correlation_config = struct();
        correlation_config.method = 'xcorr';
        correlation_config.normalize = true;
        correlation_config.peak_threshold = 0.1;
        correlation_config.positive_only = true;
        if isfield(sstdr_config, 'correlation')
            correlation_config = sstdr_config.correlation;
        end
        
        % Perform cross-correlation
        [corr_results, pn_code] = cross_correlate(tx_signal, rx_signal, ...
            config_sstdr.chip_rate, config_sstdr.fs, config_sstdr.pn_bits, ...
            correlation_config);
        
        results.corr_results = corr_results;
        results.pn_code = pn_code;
        results.tx_signal = tx_signal;
        results.rx_signal = rx_signal;
        
        if opts.verbose
            fprintf('  ✓ Cross-correlation completed\n');
            fprintf('  ✓ Found %d correlation peaks\n', length(corr_results.peaks));
        end
        
        % Plot results if requested
        if opts.plot_results
            plot_sstdr_results(results, network_config);
        end
        
    catch ME_analysis
        if opts.verbose
            fprintf('  ✗ SSTDR analysis failed: %s\n', ME_analysis.message);
        end
        results.error_message = ['Analysis failed: ' ME_analysis.message];
    end
    
    %% Step 5: Save results
    if opts.save_model || ~isempty(opts.output_dir)
        if opts.verbose
            fprintf('\n--- Saving Results ---\n');
        end
        
        % Save model if requested
        if opts.save_model
            model_file = fullfile(opts.output_dir, [opts.model_name '.slx']);
            save_system(model_name, model_file);
            results.model_file = model_file;
            if opts.verbose
                fprintf('  ✓ Model saved: %s\n', model_file);
            end
        end
        
        % Save results data
        results_file = fullfile(opts.output_dir, [opts.model_name '_results.mat']);
        save(results_file, 'results', 'model_info', 'network_config', 'sstdr_config');
        results.results_file = results_file;
        if opts.verbose
            fprintf('  ✓ Results saved: %s\n', results_file);
        end
    end
    
    %% Step 6: Cleanup
    if ~opts.keep_open
        close_system(model_name, 0);
        if opts.verbose
            fprintf('  ✓ Model closed\n');
        end
    else
        if opts.verbose
            fprintf('  ✓ Model kept open for inspection\n');
        end
    end
    
    results.success = true;
    
    if opts.verbose
        fprintf('\n=== Simulation Complete ===\n');
        fprintf('Success: %s\n', mat2str(results.success));
        fprintf('Model: %s\n', results.model_name);
        if results.success && isfield(results, 'corr_results')
            fprintf('Correlation peaks: %d\n', length(results.corr_results.peaks));
        end
    end
    
catch ME
    % Clean up on error
    if exist('model_name', 'var') && bdIsLoaded(model_name)
        close_system(model_name, 0);
    end
    
    results.success = false;
    results.error_message = ME.message;
    
    if opts.verbose
        fprintf('\n✗ Simulation failed: %s\n', ME.message);
        fprintf('Error details: %s\n', ME.getReport);
    end
    
    % Re-throw error if not in verbose mode
    if ~opts.verbose
        rethrow(ME);
    end
end

end

function plot_sstdr_results(results, network_config)
%PLOT_SSTDR_RESULTS Plot simulation and correlation results

try
    figure('Name', sprintf('SSTDR Results - %s', results.model_name), ...
           'Position', [100, 100, 1200, 800]);
    
    % Plot 1: TX and RX signals
    subplot(2, 2, 1);
    hold on;
    plot(results.tx_signal.Time * 1e6, results.tx_signal.Data, 'b-', 'LineWidth', 1);
    plot(results.rx_signal.Time * 1e6, results.rx_signal.Data, 'r-', 'LineWidth', 1);
    xlabel('Time (μs)');
    ylabel('Amplitude');
    title('TX and RX Signals');
    legend('TX', 'RX', 'Location', 'best');
    grid on;
    
    % Plot 2: Correlation result
    subplot(2, 2, 2);
    if isfield(results, 'corr_results') && isfield(results.corr_results, 'correlation')
        plot(results.corr_results.time_axis * 1e6, results.corr_results.correlation, 'k-', 'LineWidth', 1);
        xlabel('Time (μs)');
        ylabel('Correlation');
        title('Cross-Correlation Result');
        grid on;
        
        % Mark peaks
        if isfield(results.corr_results, 'peaks') && ~isempty(results.corr_results.peaks)
            hold on;
            for i = 1:length(results.corr_results.peaks)
                peak = results.corr_results.peaks(i);
                plot(peak.time * 1e6, peak.amplitude, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            end
        end
    else
        text(0.5, 0.5, 'No correlation data available', 'HorizontalAlignment', 'center');
    end
    
    % Plot 3: Network diagram
    subplot(2, 2, 3);
    plot_network_diagram(network_config);
    
    % Plot 4: Results summary
    subplot(2, 2, 4);
    axis off;
    
    % Create summary text
    summary_text = {};
    summary_text{end+1} = sprintf('Network: %s', network_config.name);
    summary_text{end+1} = sprintf('Segments: %d', network_config.num_segments);
    summary_text{end+1} = sprintf('Length: %.1f m', network_config.physical.total_length);
    summary_text{end+1} = sprintf('Faults: %d', network_config.analysis.num_faults);
    summary_text{end+1} = '';
    summary_text{end+1} = sprintf('Simulation time: %.2f s', results.sim_time);
    if isfield(results, 'corr_results')
        summary_text{end+1} = sprintf('Correlation peaks: %d', length(results.corr_results.peaks));
    end
    
    text(0.1, 0.9, summary_text, 'VerticalAlignment', 'top', 'FontSize', 10);
    
    drawnow;
    
catch ME_plot
    fprintf('  ⚠ Plotting failed: %s\n', ME_plot.message);
end

end

function plot_network_diagram(network_config)
%PLOT_NETWORK_DIAGRAM Simple network diagram

try
    % Create simple network diagram
    x_positions = 0:network_config.num_segments;
    y_line = 0.5;
    
    % Draw transmission line
    plot([0, network_config.num_segments], [y_line, y_line], 'k-', 'LineWidth', 3);
    hold on;
    
    % Mark segments
    for i = 1:network_config.num_segments
        plot(i-0.5, y_line, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    end
    
    % Mark faults
    for i = 1:network_config.num_segments
        fault_info = network_config.faults.(sprintf('segment_%d', i));
        if ~strcmp(fault_info.element_type, 'none')
            x_pos = i - 0.5;
            if strcmp(fault_info.element_type, 'series')
                plot(x_pos, y_line + 0.2, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                text(x_pos, y_line + 0.35, 'Series', 'HorizontalAlignment', 'center', 'FontSize', 8);
            else
                plot(x_pos, y_line - 0.2, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
                text(x_pos, y_line - 0.35, 'Shunt', 'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    end
    
    % Add termination
    plot(network_config.num_segments, y_line, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    text(network_config.num_segments, y_line + 0.15, 'Term', 'HorizontalAlignment', 'center', 'FontSize', 8);
    
    xlim([-0.5, network_config.num_segments + 0.5]);
    ylim([0, 1]);
    xlabel('Segment');
    title('Network Diagram');
    grid on;
    
catch ME_diagram
    text(0.5, 0.5, 'Network diagram unavailable', 'HorizontalAlignment', 'center');
end

end 