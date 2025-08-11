function [stats, fig_handles] = visualize_sstdr_dataset(dataset_path, varargin)
%VISUALIZE_SSTDR_DATASET Visualize generated SSTDR dataset with plots and statistics
%
% USAGE:
%   [stats, fig_handles] = visualize_sstdr_dataset(dataset_path)
%   [stats, fig_handles] = visualize_sstdr_dataset(dataset_path, 'network_id', 1)
%   [stats, fig_handles] = visualize_sstdr_dataset(dataset_path, 'show_all', true)
%
% INPUTS:
%   dataset_path - Path to dataset folder or specific network .mat file
%
% OPTIONAL PARAMETERS:
%   'network_id'     - Specific network to visualize (default: 1)
%   'show_all'       - Show all networks in dataset (default: false)
%   'save_plots'     - Save plots to files (default: false)
%   'plot_format'    - Format for saved plots: 'png', 'pdf', 'fig' (default: 'png')
%   'show_summary'   - Show dataset summary statistics (default: true)
%   'verbose'        - Enable verbose output (default: true)
%
% OUTPUTS:
%   stats        - Structure with dataset statistics
%   fig_handles  - Array of figure handles created
%
% EXAMPLES:
%   % Visualize network 1 from a dataset
%   [stats, figs] = visualize_sstdr_dataset('datasets/test_dataset_2025-07-17_11-51-20');
%
%   % Visualize specific network
%   [stats, figs] = visualize_sstdr_dataset('datasets/test_dataset_2025-07-17_11-51-20', ...
%       'network_id', 3, 'save_plots', true);
%
%   % Show summary of entire dataset
%   [stats, figs] = visualize_sstdr_dataset('datasets/test_dataset_2025-07-17_11-51-20', ...
%       'show_all', true, 'show_summary', true);

% Input validation
if nargin < 1
    error('Dataset path is required');
end

% Parse optional parameters
p = inputParser;
addRequired(p, 'dataset_path', @(x) ischar(x) || isstring(x));
addParameter(p, 'network_id', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'show_all', false, @islogical);
addParameter(p, 'save_plots', false, @islogical);
addParameter(p, 'plot_format', 'png', @(x) ismember(x, {'png', 'pdf', 'fig'}));
addParameter(p, 'show_summary', true, @islogical);
addParameter(p, 'verbose', true, @islogical);

parse(p, dataset_path, varargin{:});
opts = p.Results;

% Initialize outputs
stats = struct();
fig_handles = [];

% Determine if input is dataset folder or specific file
if isfile(dataset_path) && endsWith(dataset_path, '.mat')
    % Single network file
    network_data = load(dataset_path);
    dataset_folder = fileparts(dataset_path);
    single_network = true;
    if opts.verbose
        fprintf('=== Visualizing Single Network ===\n');
        fprintf('File: %s\n', dataset_path);
    end
else
    % Dataset folder
    dataset_folder = dataset_path;
    single_network = false;
    if opts.verbose
        fprintf('=== Visualizing SSTDR Dataset ===\n');
        fprintf('Dataset: %s\n', dataset_folder);
    end
end

% Get dataset statistics
if opts.show_summary || ~single_network
    stats = analyze_dataset_statistics(dataset_folder, opts.verbose);
end

% Visualize networks
if single_network
    % Visualize single network
    fig_handles = visualize_single_network(network_data.dataset, dataset_path, opts);
    
elseif opts.show_all
    % Visualize all networks in dataset
    network_files = get_network_files(dataset_folder);
    if opts.verbose
        fprintf('\nVisualizing %d networks...\n', length(network_files));
    end
    
    for i = 1:length(network_files)
        network_data = load(network_files{i});
        if opts.verbose
            fprintf('Network %d: %s\n', i, network_files{i});
        end
        figs = visualize_single_network(network_data.dataset, network_files{i}, opts);
        fig_handles = [fig_handles, figs];
    end
    
else
    % Visualize specific network
    network_files = get_network_files(dataset_folder);
    if opts.network_id > length(network_files)
        error('Network ID %d not found. Dataset has %d networks.', opts.network_id, length(network_files));
    end
    
    network_data = load(network_files{opts.network_id});
    if opts.verbose
        fprintf('\nVisualizing Network %d\n', opts.network_id);
    end
    fig_handles = visualize_single_network(network_data.dataset, network_files{opts.network_id}, opts);
end

if opts.verbose
    fprintf('\nVisualization complete! Generated %d figures.\n', length(fig_handles));
end

end

function stats = analyze_dataset_statistics(dataset_folder, verbose)
%ANALYZE_DATASET_STATISTICS Compute comprehensive dataset statistics

% Initialize statistics structure
stats = struct();
stats.dataset_folder = dataset_folder;
stats.num_networks = 0;
stats.segments = struct('min', inf, 'max', 0, 'mean', 0, 'std', 0);
stats.length = struct('min', inf, 'max', 0, 'mean', 0, 'std', 0);
stats.faults = struct('total', 0, 'per_network_mean', 0, 'per_network_std', 0);
stats.fault_types = struct('series', 0, 'shunt', 0);
stats.termination = struct('min', inf, 'max', 0, 'mean', 0, 'std', 0);

% Get all network files
network_files = get_network_files(dataset_folder);
stats.num_networks = length(network_files);

if stats.num_networks == 0
    if verbose
        fprintf('No networks found in dataset folder.\n');
    end
    return;
end

% Collect data from all networks
segments_data = zeros(stats.num_networks, 1);
length_data = zeros(stats.num_networks, 1);
faults_data = zeros(stats.num_networks, 1);
termination_data = zeros(stats.num_networks, 1);

for i = 1:stats.num_networks
    try
        network_data = load(network_files{i});
        dataset = network_data.dataset;
        
        % Network configuration
        if isfield(dataset, 'network_config')
            nc = dataset.network_config;
            
            % Segments
            if isfield(nc, 'num_segments')
                segments_data(i) = nc.num_segments;
            elseif isfield(nc, 'load_vector')
                segments_data(i) = length(nc.load_vector);
            end
            
            % Total length
            if isfield(nc, 'physical') && isfield(nc.physical, 'total_length')
                length_data(i) = nc.physical.total_length;
            elseif isfield(nc, 'dx') && segments_data(i) > 0
                length_data(i) = nc.dx * segments_data(i);
            end
            
            % Count faults
            if isfield(nc, 'load_vector')
                load_vec = nc.load_vector;
                fault_count = sum(load_vec ~= 0);
                faults_data(i) = fault_count;
                
                % Count fault types (positive = series, negative = shunt)
                stats.fault_types.series = stats.fault_types.series + sum(load_vec > 0);
                stats.fault_types.shunt = stats.fault_types.shunt + sum(load_vec < 0);
            end
            
            % Termination impedance
            if isfield(nc, 'physical') && isfield(nc.physical, 'termination_impedance')
                termination_data(i) = nc.physical.termination_impedance;
            end
        end
        
    catch ME
        if verbose
            fprintf('Warning: Could not load network %d: %s\n', i, ME.message);
        end
    end
end

% Compute statistics
valid_segments = segments_data(segments_data > 0);
valid_lengths = length_data(length_data > 0);
valid_terminations = termination_data(termination_data > 0);

if ~isempty(valid_segments)
    stats.segments.min = min(valid_segments);
    stats.segments.max = max(valid_segments);
    stats.segments.mean = mean(valid_segments);
    stats.segments.std = std(valid_segments);
end

if ~isempty(valid_lengths)
    stats.length.min = min(valid_lengths);
    stats.length.max = max(valid_lengths);
    stats.length.mean = mean(valid_lengths);
    stats.length.std = std(valid_lengths);
end

if ~isempty(valid_terminations)
    stats.termination.min = min(valid_terminations);
    stats.termination.max = max(valid_terminations);
    stats.termination.mean = mean(valid_terminations);
    stats.termination.std = std(valid_terminations);
end

stats.faults.total = sum(faults_data);
stats.faults.per_network_mean = mean(faults_data);
stats.faults.per_network_std = std(faults_data);

% Display summary
if verbose
    fprintf('\n=== Dataset Statistics ===\n');
    fprintf('Networks: %d\n', stats.num_networks);
    fprintf('Segments per network: %.1f ± %.1f (range: %d - %d)\n', ...
        stats.segments.mean, stats.segments.std, stats.segments.min, stats.segments.max);
    fprintf('Network length: %.1f ± %.1f m (range: %.1f - %.1f m)\n', ...
        stats.length.mean, stats.length.std, stats.length.min, stats.length.max);
    fprintf('Faults per network: %.1f ± %.1f (total: %d)\n', ...
        stats.faults.per_network_mean, stats.faults.per_network_std, stats.faults.total);
    fprintf('Fault types: %d series, %d shunt\n', ...
        stats.fault_types.series, stats.fault_types.shunt);
    fprintf('Termination impedance: %.1f ± %.1f Ω (range: %.1f - %.1f Ω)\n', ...
        stats.termination.mean, stats.termination.std, stats.termination.min, stats.termination.max);
end

end

function fig_handles = visualize_single_network(dataset, file_path, opts)
%VISUALIZE_SINGLE_NETWORK Create visualizations for a single network

fig_handles = [];

% Extract network information
if isfield(dataset, 'network_config')
    network_config = dataset.network_config;
else
    error('Network configuration not found in dataset');
end

% Get network name for titles
if isfield(dataset, 'metadata') && isfield(dataset.metadata, 'model_name')
    network_name = dataset.metadata.model_name;
else
    [~, network_name, ~] = fileparts(file_path);
end

% Create main visualization figure
fig1 = figure('Name', sprintf('SSTDR Analysis: %s', network_name), ...
    'Position', [100, 100, 1200, 800]);
fig_handles(end+1) = fig1;

% Create subplots
subplot(2, 3, [1, 2]); % Network diagram (spans 2 columns)
plot_network_diagram(network_config);
title('Network Topology and Faults');

subplot(2, 3, 3); % Network statistics
plot_network_stats(dataset);
title('Network Statistics');

subplot(2, 3, 4); % TX/RX signals
plot_signals(dataset);
title('SSTDR Signals');

subplot(2, 3, 5); % Correlation results
plot_correlation(dataset);
title('Cross-Correlation');

subplot(2, 3, 6); % Fault analysis
plot_fault_analysis(dataset);
title('Fault Analysis');

% Add overall title
sgtitle(sprintf('SSTDR Dataset Visualization: %s', network_name), 'FontSize', 14, 'FontWeight', 'bold');

% Save plots if requested
if opts.save_plots
    save_plot_file(fig1, file_path, 'overview', opts.plot_format);
end

% Create detailed correlation figure if correlation data exists
if isfield(dataset, 'correlation_results')
    fig2 = figure('Name', sprintf('Correlation Details: %s', network_name), ...
        'Position', [150, 150, 1000, 600]);
    fig_handles(end+1) = fig2;
    
    plot_detailed_correlation(dataset);
    sgtitle(sprintf('Detailed Cross-Correlation Analysis: %s', network_name), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    if opts.save_plots
        save_plot_file(fig2, file_path, 'correlation', opts.plot_format);
    end
end

end

function plot_network_diagram(network_config)
%PLOT_NETWORK_DIAGRAM Visualize network topology and faults

if ~isfield(network_config, 'load_vector')
    text(0.5, 0.5, 'Network diagram not available', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

load_vector = network_config.load_vector;
num_segments = length(load_vector);

% Get segment length
if isfield(network_config, 'dx')
    dx = network_config.dx;
elseif isfield(network_config, 'physical') && isfield(network_config.physical, 'total_length')
    dx = network_config.physical.total_length / num_segments;
else
    dx = 10; % Default segment length
end

% Create position array
positions = (0:num_segments) * dx;

% Plot transmission line
line([0, positions(end)], [0, 0], 'Color', 'k', 'LineWidth', 3);
hold on;

% Plot segments
for i = 1:num_segments
    % Segment markers
    plot(positions(i), 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    
    % Segment labels
    text(positions(i) + dx/2, 0.15, sprintf('S%d', i), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
    
    % Fault visualization
    fault_val = load_vector(i);
    if fault_val ~= 0
        if fault_val > 0
            % Series fault (resistance symbol)
            plot([positions(i) + dx*0.3, positions(i) + dx*0.7], [0, 0], 'r-', 'LineWidth', 2);
            plot([positions(i) + dx*0.4, positions(i) + dx*0.6], [0.1, -0.1], 'r-', 'LineWidth', 2);
            plot([positions(i) + dx*0.4, positions(i) + dx*0.6], [-0.1, 0.1], 'r-', 'LineWidth', 2);
            text(positions(i) + dx/2, -0.25, sprintf('R=%.0fΩ', abs(fault_val)*1000), ...
                'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'r');
        else
            % Shunt fault (capacitor/short symbol)
            plot(positions(i) + dx/2, 0, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            plot([positions(i) + dx/2, positions(i) + dx/2], [0, -0.2], 'r-', 'LineWidth', 2);
            text(positions(i) + dx/2, -0.35, sprintf('G=%.0fΩ', abs(fault_val)*1000), ...
                'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'r');
        end
    end
end

% Plot termination
plot(positions(end), 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
if isfield(network_config, 'physical') && isfield(network_config.physical, 'termination_impedance')
    term_z = network_config.physical.termination_impedance;
    text(positions(end), -0.15, sprintf('Z_L=%.1fΩ', term_z), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

% Add source
plot(0, 0, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
text(0, 0.15, 'SSTDR', 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');

% Formatting
xlim([-dx*0.2, positions(end) + dx*0.2]);
ylim([-0.5, 0.3]);
xlabel('Distance (m)');
ylabel('');
grid on;
axis equal;
hold off;

% Add legend
legend_elements = {'Transmission Line', 'Series Fault', 'Shunt Fault', 'Source', 'Termination'};
legend(legend_elements, 'Location', 'northeast', 'FontSize', 8);

end

function plot_network_stats(dataset)
%PLOT_NETWORK_STATS Display network statistics as text

axis off;

stats_text = {};
if isfield(dataset, 'network_config')
    nc = dataset.network_config;
    
    % Basic network info
    if isfield(nc, 'load_vector')
        load_vec = nc.load_vector;
        num_segments = length(load_vec);
        num_faults = sum(load_vec ~= 0);
        num_series = sum(load_vec > 0);
        num_shunt = sum(load_vec < 0);
        
        stats_text{end+1} = sprintf('Segments: %d', num_segments);
        stats_text{end+1} = sprintf('Total Faults: %d', num_faults);
        stats_text{end+1} = sprintf('Series Faults: %d', num_series);
        stats_text{end+1} = sprintf('Shunt Faults: %d', num_shunt);
    end
    
    % Physical parameters
    if isfield(nc, 'physical')
        phys = nc.physical;
        if isfield(phys, 'total_length')
            stats_text{end+1} = sprintf('Length: %.1f m', phys.total_length);
        end
        if isfield(phys, 'Z0')
            stats_text{end+1} = sprintf('Z₀: %.1f Ω', phys.Z0);
        end
        if isfield(phys, 'termination_impedance')
            stats_text{end+1} = sprintf('Z_L: %.1f Ω', phys.termination_impedance);
        end
    end
end

% Simulation info
if isfield(dataset, 'simulation_config')
    sc = dataset.simulation_config;
    if isfield(sc, 'pn_config')
        pn = sc.pn_config;
        if isfield(pn, 'chip_rate_hz')
            stats_text{end+1} = sprintf('Chip Rate: %.0f kHz', pn.chip_rate_hz/1000);
        end
        if isfield(pn, 'carrier_freq_hz')
            stats_text{end+1} = sprintf('Carrier: %.0f kHz', pn.carrier_freq_hz/1000);
        end
    end
end

% Performance info
if isfield(dataset, 'performance')
    perf = dataset.performance;
    if isfield(perf, 'simulation_time_s')
        stats_text{end+1} = sprintf('Sim Time: %.2f s', perf.simulation_time_s);
    end
end

% Display text
for i = 1:length(stats_text)
    text(0.05, 0.95 - (i-1)*0.08, stats_text{i}, 'FontSize', 10, 'FontWeight', 'normal');
end

xlim([0, 1]);
ylim([0, 1]);

end

function plot_signals(dataset)
%PLOT_SIGNALS Plot TX and RX signals

if isfield(dataset, 'tx_signal') && isfield(dataset, 'rx_signal')
    tx = dataset.tx_signal;
    rx = dataset.rx_signal;
    
    % Plot TX signal
    if isfield(tx, 'Time') && isfield(tx, 'Data')
        plot(tx.Time * 1000, tx.Data, 'b-', 'LineWidth', 1);
        hold on;
    end
    
    % Plot RX signal
    if isfield(rx, 'Time') && isfield(rx, 'Data')
        plot(rx.Time * 1000, rx.Data, 'r-', 'LineWidth', 1);
    end
    
    xlabel('Time (ms)');
    ylabel('Amplitude (V)');
    legend('TX Signal', 'RX Signal', 'Location', 'best');
    grid on;
    hold off;
else
    text(0.5, 0.5, 'Signal data not available', 'HorizontalAlignment', 'center');
    axis off;
end

end

function plot_correlation(dataset)
%PLOT_CORRELATION Plot correlation results

if isfield(dataset, 'correlation_results')
    corr = dataset.correlation_results;
    
    if isfield(corr, 'correlation') && isfield(corr, 'lags')
        % Convert lags to time (assuming sampling rate info is available)
        if isfield(corr, 'config') && isfield(corr.config, 'fs')
            fs = corr.config.fs;
            time_lags = corr.lags / fs * 1000; % Convert to ms
        else
            time_lags = corr.lags;
        end
        
        plot(time_lags, abs(corr.correlation), 'k-', 'LineWidth', 1);
        xlabel('Time Lag (ms)');
        ylabel('|Correlation|');
        grid on;
        
        % Mark peaks if available
        if isfield(corr, 'peaks')
            hold on;
            peaks = corr.peaks;
            if isfield(peaks, 'indices') && isfield(peaks, 'values')
                peak_times = time_lags(peaks.indices);
                plot(peak_times, abs(peaks.values), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
            end
            hold off;
        end
    else
        text(0.5, 0.5, 'Correlation data not available', 'HorizontalAlignment', 'center');
        axis off;
    end
else
    text(0.5, 0.5, 'Correlation results not available', 'HorizontalAlignment', 'center');
    axis off;
end

end

function plot_fault_analysis(dataset)
%PLOT_FAULT_ANALYSIS Analyze and display fault-related information

axis off;

fault_text = {};

if isfield(dataset, 'network_config') && isfield(dataset.network_config, 'load_vector')
    load_vec = dataset.network_config.load_vector;
    
    % Find fault locations
    fault_indices = find(load_vec ~= 0);
    
    if ~isempty(fault_indices)
        fault_text{end+1} = 'Fault Locations:';
        
        for i = 1:length(fault_indices)
            idx = fault_indices(i);
            fault_val = load_vec(idx);
            
            if fault_val > 0
                fault_type = 'Series';
                resistance = abs(fault_val) * 1000; % Convert to Ohms
            else
                fault_type = 'Shunt';
                resistance = abs(fault_val) * 1000; % Convert to Ohms
            end
            
            fault_text{end+1} = sprintf('Seg %d: %s %.0fΩ', idx, fault_type, resistance);
        end
    else
        fault_text{end+1} = 'No faults detected';
    end
    
    % Calculate fault density
    fault_density = length(fault_indices) / length(load_vec) * 100;
    fault_text{end+1} = '';
    fault_text{end+1} = sprintf('Fault Density: %.1f%%', fault_density);
    
else
    fault_text{end+1} = 'Fault analysis not available';
end

% Display text
for i = 1:length(fault_text)
    text(0.05, 0.95 - (i-1)*0.1, fault_text{i}, 'FontSize', 9, 'FontWeight', 'normal');
end

xlim([0, 1]);
ylim([0, 1]);

end

function plot_detailed_correlation(dataset)
%PLOT_DETAILED_CORRELATION Create detailed correlation analysis plots

if ~isfield(dataset, 'correlation_results')
    return;
end

corr = dataset.correlation_results;

% Subplot 1: Full correlation
subplot(2, 2, 1);
if isfield(corr, 'correlation') && isfield(corr, 'lags')
    plot(corr.lags, real(corr.correlation), 'b-', 'LineWidth', 1);
    hold on;
    plot(corr.lags, imag(corr.correlation), 'r-', 'LineWidth', 1);
    legend('Real', 'Imaginary');
    xlabel('Lag Samples');
    ylabel('Correlation');
    title('Full Cross-Correlation');
    grid on;
    hold off;
end

% Subplot 2: Magnitude spectrum
subplot(2, 2, 2);
if isfield(corr, 'correlation')
    plot(corr.lags, abs(corr.correlation), 'k-', 'LineWidth', 1);
    xlabel('Lag Samples');
    ylabel('|Correlation|');
    title('Correlation Magnitude');
    grid on;
    
    % Mark peaks
    if isfield(corr, 'peaks') && isfield(corr.peaks, 'indices')
        hold on;
        peak_indices = corr.peaks.indices;
        peak_values = abs(corr.correlation(peak_indices));
        plot(corr.lags(peak_indices), peak_values, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        hold off;
    end
end

% Subplot 3: Peak analysis
subplot(2, 2, 3);
if isfield(corr, 'peaks')
    peaks = corr.peaks;
    if isfield(peaks, 'indices') && isfield(peaks, 'values')
        stem(peaks.indices, abs(peaks.values), 'r-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel('Peak Index');
        ylabel('Peak Magnitude');
        title('Detected Peaks');
        grid on;
    end
end

% Subplot 4: Configuration info
subplot(2, 2, 4);
axis off;
config_text = {};

if isfield(corr, 'config')
    cfg = corr.config;
    config_text{end+1} = 'Correlation Config:';
    
    if isfield(cfg, 'method')
        config_text{end+1} = sprintf('Method: %s', cfg.method);
    end
    if isfield(cfg, 'peak_threshold')
        config_text{end+1} = sprintf('Peak Threshold: %.3f', cfg.peak_threshold);
    end
    if isfield(cfg, 'normalize')
        config_text{end+1} = sprintf('Normalized: %s', string(cfg.normalize));
    end
    if isfield(cfg, 'positive_only')
        config_text{end+1} = sprintf('Positive Only: %s', string(cfg.positive_only));
    end
end

% Display configuration text
for i = 1:length(config_text)
    text(0.05, 0.95 - (i-1)*0.12, config_text{i}, 'FontSize', 10);
end

end

function network_files = get_network_files(dataset_folder)
%GET_NETWORK_FILES Find all network .mat files in dataset folder

network_files = {};

% Look for network subdirectories
network_dirs = dir(fullfile(dataset_folder, 'network_*'));
network_dirs = network_dirs([network_dirs.isdir]);

for i = 1:length(network_dirs)
    network_dir = fullfile(dataset_folder, network_dirs(i).name);
    mat_files = dir(fullfile(network_dir, '*.mat'));
    
    for j = 1:length(mat_files)
        network_files{end+1} = fullfile(network_dir, mat_files(j).name);
    end
end

% Sort network files
network_files = sort(network_files);

end

function save_plot_file(fig_handle, original_path, suffix, format)
%SAVE_PLOT_FILE Save figure to file with appropriate naming

[path_dir, name, ~] = fileparts(original_path);
plot_filename = sprintf('%s_%s.%s', name, suffix, format);
plot_path = fullfile(path_dir, plot_filename);

switch format
    case 'png'
        saveas(fig_handle, plot_path, 'png');
    case 'pdf'
        saveas(fig_handle, plot_path, 'pdf');
    case 'fig'
        savefig(fig_handle, plot_path);
end

fprintf('Saved plot: %s\n', plot_path);

end 