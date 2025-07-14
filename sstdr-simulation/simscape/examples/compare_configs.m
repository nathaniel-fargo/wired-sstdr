%% COMPARE SSTDR CONFIGURATIONS
% Runs multiple configurations and compares performance

clear; clc; close all;

%% Setup paths for organized folder structure
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'functions', 'network_generation'));
addpath(fullfile(parent_dir, 'functions', 'sstdr_simulation'));
addpath(fullfile(parent_dir, 'config'));

%% 🎯 CHANGE THIS TO YOUR MODEL NAME
MODEL_NAME = '';  % Set to your Simscape model name, or leave empty

%% Configuration comparison setup
configs_to_test = {
    'default',      % 100 kHz carrier, 1 MHz sampling
    'sine_fast',    % 250 kHz carrier, 2 MHz sampling
    'high_res',     % 100 kHz carrier, 5 MHz sampling
    'unmodulated'   % No carrier, 1 MHz sampling
};

% Storage for results
comparison_results = cell(length(configs_to_test), 1);
config_names = cell(length(configs_to_test), 1);

fprintf('=== SSTDR Configuration Comparison ===\n');
fprintf('Testing %d configurations...\n\n', length(configs_to_test));

%% Run each configuration
for i = 1:length(configs_to_test)
    config_name = configs_to_test{i};
    
    fprintf('--- Configuration %d/%d: %s ---\n', i, length(configs_to_test), config_name);
    
    % Load configuration
    sstdr_config(config_name);
    config = evalin('base', 'sstdr_config');
    config_names{i} = config.name;
    
    % Configure model if specified
    if ~isempty(MODEL_NAME)
        try
            configure_model_sampling(MODEL_NAME);
            run_sstdr_analysis('skip', MODEL_NAME, 'method', 'freq', 'plot_results', false);
        catch ME
            fprintf('⚠ Analysis failed: %s\n', ME.message);
            continue;
        end
    else
        % Use existing data
        try
            cross_correlate('method', 'freq', 'plot_results', false);
        catch ME
            fprintf('⚠ Analysis failed: %s\n', ME.message);
            continue;
        end
    end
    
    % Store results
    if evalin('base', 'exist(''correlation_results'', ''var'')')
        comparison_results{i} = evalin('base', 'correlation_results');
        fprintf('✓ Configuration completed\n');
    else
        fprintf('✗ No results generated\n');
    end
    
    fprintf('\n');
end

%% Display comparison summary
fprintf('=== COMPARISON SUMMARY ===\n');
fprintf('%-20s %-15s %-15s %-15s %-10s\n', 'Configuration', 'Modulation', 'Carrier (kHz)', 'Sampling (kHz)', 'Time (ms)');
fprintf('%s\n', repmat('-', 1, 80));

for i = 1:length(comparison_results)
    if ~isempty(comparison_results{i})
        result = comparison_results{i};
        config = result.config;
        
        % Extract configuration info
        modulation = config.pn_config.modulation;
        if strcmp(modulation, 'none')
            carrier_str = 'N/A';
        else
            carrier_str = sprintf('%.0f', config.pn_config.carrier_freq/1000);
        end
        sampling_str = sprintf('%.0f', config.pn_config.fs/1000);
        
        % Get computation time
        if isfield(result, 'compute_time_freq')
            time_str = sprintf('%.1f', result.compute_time_freq*1000);
        else
            time_str = 'N/A';
        end
        
        fprintf('%-20s %-15s %-15s %-15s %-10s\n', ...
                config_names{i}, modulation, carrier_str, sampling_str, time_str);
        
        % Show detected peaks
        if isfield(result, 'peak_times_freq') && ~isempty(result.peak_times_freq)
            fprintf('  Peaks at: ');
            fprintf('%.1f ', result.peak_times_freq * 1e6);
            fprintf('µs\n');
        end
    else
        fprintf('%-20s %-15s %-15s %-15s %-10s\n', config_names{i}, 'FAILED', '-', '-', '-');
    end
end

%% Create comparison plots
if sum(~cellfun(@isempty, comparison_results)) >= 2
    fprintf('\nCreating comparison plots...\n');
    
    figure('Position', [100, 100, 1400, 1000]);
    
    valid_results = comparison_results(~cellfun(@isempty, comparison_results));
    valid_names = config_names(~cellfun(@isempty, comparison_results));
    
    num_configs = length(valid_results);
    
    for i = 1:num_configs
        subplot(2, ceil(num_configs/2), i);
        
        result = valid_results{i};
        if isfield(result, 'correlation_freq') && isfield(result, 'time_lags_freq')
            plot(result.time_lags_freq * 1e6, real(result.correlation_freq), 'LineWidth', 2);
            title(sprintf('%s', valid_names{i}), 'Interpreter', 'none');
            xlabel('Time Lag (µs)');
            ylabel('Correlation');
            grid on;
            
            % Mark peaks
            if isfield(result, 'peak_times_freq') && ~isempty(result.peak_times_freq)
                hold on;
                peak_indices = result.peak_locs_freq;
                plot(result.peak_times_freq * 1e6, real(result.correlation_freq(peak_indices)), ...
                     'ro', 'MarkerSize', 8, 'LineWidth', 2);
                hold off;
            end
        end
    end
    
    sgtitle('SSTDR Configuration Comparison', 'FontSize', 16, 'FontWeight', 'bold');
end

%% Performance ranking
fprintf('\n=== PERFORMANCE RANKING ===\n');
valid_times = [];
valid_indices = [];

for i = 1:length(comparison_results)
    if ~isempty(comparison_results{i}) && isfield(comparison_results{i}, 'compute_time_freq')
        valid_times(end+1) = comparison_results{i}.compute_time_freq * 1000; %#ok<AGROW>
        valid_indices(end+1) = i; %#ok<AGROW>
    end
end

if ~isempty(valid_times)
    [sorted_times, sort_idx] = sort(valid_times);
    
    fprintf('Fastest to slowest correlation:\n');
    for i = 1:length(sorted_times)
        original_idx = valid_indices(sort_idx(i));
        fprintf('%d. %s: %.1f ms\n', i, config_names{original_idx}, sorted_times(i));
    end
end

fprintf('\n🎉 Comparison complete! Check the plots for visual comparison.\n'); 