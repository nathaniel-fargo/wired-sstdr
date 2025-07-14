function analyze_sstdr_dataset(dataset_path, varargin)
%ANALYZE_SSTDR_DATASET Analyze and visualize SSTDR dataset results
%
% Usage:
%   analyze_sstdr_dataset('datasets/my_dataset')              % Analyze all networks
%   analyze_sstdr_dataset('datasets/my_dataset', 'display', 'single')  % Single display window
%   analyze_sstdr_dataset('datasets/my_dataset', 'display', 'multiple')  % Multiple windows
%   analyze_sstdr_dataset('datasets/my_dataset', 'display', 'first_n', 'n', 5)  % First 5 networks
%   analyze_sstdr_dataset('datasets/my_dataset', 'save_plots', true)  % Save plots to files
%
% Parameters:
%   dataset_path - Path to dataset directory or summary file
%   'display' - Display mode: 'single', 'multiple', 'first_n', 'none' (default: 'single')
%   'n' - Number of networks to display when display='first_n' (default: 5)
%   'save_plots' - Save plots to files (default: false)
%   'plot_types' - Cell array of plot types: {'correlation', 'waveform', 'network'} (default: all)
%   'verbose' - Display analysis information (default: true)

%% Parse inputs
if nargin < 1
    error('Dataset path is required');
end

p = inputParser;
addRequired(p, 'dataset_path', @ischar);
addParameter(p, 'display', 'single', @(x) ismember(x, {'single', 'multiple', 'first_n', 'none'}));
addParameter(p, 'n', 5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'save_plots', false, @islogical);
addParameter(p, 'plot_types', {'correlation', 'waveform', 'network'}, @iscell);
addParameter(p, 'verbose', true, @islogical);

parse(p, dataset_path, varargin{:});
opts = p.Results;

%% Load dataset
if opts.verbose
    fprintf('=== SSTDR Dataset Analysis ===\n');
    fprintf('Loading dataset: %s\n', dataset_path);
end

% Check if path is a directory or file
if isfolder(dataset_path)
    summary_file = fullfile(dataset_path, 'dataset_summary.mat');
    if ~exist(summary_file, 'file')
        error('Dataset summary file not found: %s', summary_file);
    end
else
    summary_file = dataset_path;
    if ~exist(summary_file, 'file')
        error('Dataset file not found: %s', summary_file);
    end
end

% Load dataset
load(summary_file, 'dataset_results', 'config');

if opts.verbose
    fprintf('Dataset loaded successfully!\n');
    fprintf('Configuration: %s\n', config.metadata.name);
    fprintf('Total networks: %d\n', dataset_results.metadata.total_networks);
    fprintf('Successful networks: %d (%.1f%%)\n', ...
        dataset_results.metadata.successful_networks, ...
        dataset_results.metadata.success_rate * 100);
end

%% Determine which networks to analyze
total_networks = dataset_results.metadata.total_networks;
successful_networks = dataset_results.metadata.successful_networks;

switch opts.display
    case 'single'
        networks_to_analyze = 1:total_networks;
        num_figures = 1;
    case 'multiple'
        networks_to_analyze = 1:total_networks;
        num_figures = successful_networks;
    case 'first_n'
        networks_to_analyze = 1:min(opts.n, total_networks);
        num_figures = min(opts.n, successful_networks);
    case 'none'
        networks_to_analyze = [];
        num_figures = 0;
end

%% Analyze dataset statistics
if opts.verbose
    fprintf('\n--- Dataset Statistics ---\n');
    analyze_dataset_statistics(dataset_results, config);
end

%% Create visualizations
if num_figures > 0
    if opts.verbose
        fprintf('\n--- Creating Visualizations ---\n');
        fprintf('Display mode: %s\n', opts.display);
        fprintf('Networks to analyze: %d\n', length(networks_to_analyze));
    end
    
    create_visualizations(dataset_results, config, networks_to_analyze, opts);
end

if opts.verbose
    fprintf('\nAnalysis complete!\n');
end

end

function analyze_dataset_statistics(dataset_results, config)
%ANALYZE_DATASET_STATISTICS Display comprehensive dataset statistics

fprintf('Processing Times:\n');
fprintf('  Total time: %.1f minutes\n', dataset_results.metadata.total_time / 60);
fprintf('  Average per network: %.2f seconds\n', dataset_results.metadata.avg_time_per_network);
fprintf('  Min time: %.2f seconds\n', min(dataset_results.metadata.processing_times));
fprintf('  Max time: %.2f seconds\n', max(dataset_results.metadata.processing_times));

% Analyze successful networks
successful_count = 0;
fault_counts = [];
segment_counts = [];
correlation_peaks = [];

for i = 1:dataset_results.metadata.total_networks
    network = dataset_results.networks{i};
    
    if isfield(network, 'corr_results')
        successful_count = successful_count + 1;
        
        % Count faults
        if isfield(network.net_config, 'faults')
            fault_counts(end+1) = length(network.net_config.faults);
        end
        
        % Count segments
        if isfield(network.net_config, 'segments')
            segment_counts(end+1) = length(network.net_config.segments);
        end
        
        % Get correlation peaks
        if isfield(network.corr_results, 'peak_magnitudes')
            correlation_peaks = [correlation_peaks; network.corr_results.peak_magnitudes(:)];
        end
    end
end

fprintf('\nNetwork Statistics:\n');
if ~isempty(fault_counts)
    fprintf('  Faults per network: %.1f ± %.1f (range: %d-%d)\n', ...
        mean(fault_counts), std(fault_counts), min(fault_counts), max(fault_counts));
end

if ~isempty(segment_counts)
    fprintf('  Segments per network: %.1f ± %.1f (range: %d-%d)\n', ...
        mean(segment_counts), std(segment_counts), min(segment_counts), max(segment_counts));
end

if ~isempty(correlation_peaks)
    fprintf('  Correlation peaks: %.3f ± %.3f (range: %.3f-%.3f)\n', ...
        mean(correlation_peaks), std(correlation_peaks), min(correlation_peaks), max(correlation_peaks));
end

fprintf('\nConfiguration Summary:\n');
fprintf('  SSTDR: %.1f MHz chip rate, %.1f MHz sampling\n', config.sstdr.chip_rate/1e6, config.sstdr.fs/1e6);
fprintf('  PN Sequence: %d bits (%d chips)\n', config.sstdr.pn_bits, 2^config.sstdr.pn_bits - 1);
fprintf('  Simulation: %.1f μs duration\n', config.simulation.duration * 1e6);
fprintf('  Network: %.1f km segments, %.1f%% fault probability\n', ...
    config.network.dx/1000, config.network.fault_probability * 100);

end

function create_visualizations(dataset_results, config, networks_to_analyze, opts)
%CREATE_VISUALIZATIONS Create plots based on display mode and options

successful_networks = [];
for i = networks_to_analyze
    if i <= length(dataset_results.networks) && isfield(dataset_results.networks{i}, 'corr_results')
        successful_networks(end+1) = i;
    end
end

if isempty(successful_networks)
    fprintf('No successful networks to visualize\n');
    return;
end

switch opts.display
    case 'single'
        create_single_window_plots(dataset_results, config, successful_networks, opts);
    case 'multiple'
        create_multiple_window_plots(dataset_results, config, successful_networks, opts);
    case 'first_n'
        create_multiple_window_plots(dataset_results, config, successful_networks(1:min(opts.n, length(successful_networks))), opts);
end

end

function create_single_window_plots(dataset_results, config, networks, opts)
%CREATE_SINGLE_WINDOW_PLOTS Create all plots in a single window

num_networks = length(networks);
num_plot_types = length(opts.plot_types);

% Calculate subplot arrangement
cols = ceil(sqrt(num_networks * num_plot_types));
rows = ceil((num_networks * num_plot_types) / cols);

figure('Name', 'SSTDR Dataset Analysis - Single Window', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 100, 1200, 800]);

plot_idx = 1;

for i = 1:num_networks
    network_idx = networks(i);
    network = dataset_results.networks{network_idx};
    
    for j = 1:num_plot_types
        subplot(rows, cols, plot_idx);
        
        switch opts.plot_types{j}
            case 'correlation'
                plot_correlation(network, network_idx, config);
            case 'waveform'
                plot_waveform(network, network_idx, config);
            case 'network'
                plot_network_diagram(network, network_idx, config);
        end
        
        plot_idx = plot_idx + 1;
    end
end

if opts.save_plots
    save_path = fullfile(config.dataset.output_dir, 'analysis_single_window.png');
    saveas(gcf, save_path);
    fprintf('Saved single window plot: %s\n', save_path);
end

end

function create_multiple_window_plots(dataset_results, config, networks, opts)
%CREATE_MULTIPLE_WINDOW_PLOTS Create separate window for each network

for i = 1:length(networks)
    network_idx = networks(i);
    network = dataset_results.networks{network_idx};
    
    figure('Name', sprintf('SSTDR Network %d Analysis', network_idx), 'NumberTitle', 'off');
    set(gcf, 'Position', [100 + (i-1)*50, 100 + (i-1)*50, 1000, 600]);
    
    num_plots = length(opts.plot_types);
    
    for j = 1:num_plots
        subplot(1, num_plots, j);
        
        switch opts.plot_types{j}
            case 'correlation'
                plot_correlation(network, network_idx, config);
            case 'waveform'
                plot_waveform(network, network_idx, config);
            case 'network'
                plot_network_diagram(network, network_idx, config);
        end
    end
    
    if opts.save_plots
        save_path = fullfile(config.dataset.output_dir, sprintf('analysis_network_%03d.png', network_idx));
        saveas(gcf, save_path);
        fprintf('Saved network %d plot: %s\n', network_idx, save_path);
    end
end

end

function plot_correlation(network, network_idx, config)
%PLOT_CORRELATION Plot correlation results

if isfield(network.corr_results, 'correlation_time') && isfield(network.corr_results, 'time_lags')
    plot(network.corr_results.time_lags * 1e6, network.corr_results.correlation_time);
    xlabel('Time (μs)');
    ylabel('Correlation');
    title(sprintf('Network %d - Correlation', network_idx));
    grid on;
    
    % Mark peaks if available
    if isfield(network.corr_results, 'peak_times_time')
        hold on;
        peak_times = network.corr_results.peak_times_time * 1e6;
        peak_mags = network.corr_results.peak_magnitudes;
        plot(peak_times, peak_mags, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
        hold off;
    end
else
    text(0.5, 0.5, 'No correlation data', 'HorizontalAlignment', 'center');
    title(sprintf('Network %d - No Correlation Data', network_idx));
end

end

function plot_waveform(network, network_idx, config)
%PLOT_WAVEFORM Plot transmitted and received waveforms

if isfield(network.sim_out, 'tx_signal') && isfield(network.sim_out, 'rx_signal')
    time = network.sim_out.tx_signal.Time * 1e6;
    
    plot(time, network.sim_out.tx_signal.Data, 'b-', 'LineWidth', 1);
    hold on;
    plot(time, network.sim_out.rx_signal.Data, 'r-', 'LineWidth', 1);
    hold off;
    
    xlabel('Time (μs)');
    ylabel('Amplitude');
    title(sprintf('Network %d - Waveforms', network_idx));
    legend('TX', 'RX', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No waveform data', 'HorizontalAlignment', 'center');
    title(sprintf('Network %d - No Waveform Data', network_idx));
end

end

function plot_network_diagram(network, network_idx, config)
%PLOT_NETWORK_DIAGRAM Plot network topology diagram

if isfield(network.net_config, 'segments')
    segments = network.net_config.segments;
    num_segments = length(segments);
    
    % Create position array
    positions = 0:config.network.dx/1000:(num_segments-1)*config.network.dx/1000;
    
    % Plot transmission line
    plot(positions, zeros(size(positions)), 'k-', 'LineWidth', 3);
    hold on;
    
    % Plot faults if available
    if isfield(network.net_config, 'faults')
        faults = network.net_config.faults;
        for i = 1:length(faults)
            fault_pos = faults(i).position * config.network.dx/1000;
            fault_type = faults(i).type;
            
            if strcmp(fault_type, 'series')
                plot(fault_pos, 0, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
            else
                plot(fault_pos, 0, 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'blue');
            end
        end
    end
    
    % Plot segment boundaries
    for i = 1:num_segments
        plot([positions(i), positions(i)], [-0.1, 0.1], 'k-', 'LineWidth', 1);
    end
    
    hold off;
    
    xlabel('Distance (km)');
    ylabel('');
    title(sprintf('Network %d - Topology (%d segments)', network_idx, num_segments));
    ylim([-0.2, 0.2]);
    grid on;
    
    % Add legend
    legend_items = {'Transmission Line'};
    if isfield(network.net_config, 'faults') && ~isempty(network.net_config.faults)
        legend_items{end+1} = 'Series Fault';
        legend_items{end+1} = 'Shunt Fault';
    end
    legend(legend_items, 'Location', 'best');
    
else
    text(0.5, 0.5, 'No network data', 'HorizontalAlignment', 'center');
    title(sprintf('Network %d - No Network Data', network_idx));
end

end 