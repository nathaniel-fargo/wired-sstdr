%% Test: How Sampling Rate Affects SSTDR Signal Time Scale
% This demonstrates why changing fs dramatically changes your graphs

clear; clc; close all;

% Setup paths
current_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(current_dir, '..', 'functions', 'network_generation'));
addpath(fullfile(current_dir, '..', 'functions', 'sstdr_simulation'));
addpath(fullfile(current_dir, '..', 'config'));

%% Test different sampling rates with same carrier frequency
carrier_freq = 10e6;  % Fixed 10 MHz carrier
fs_values = [1e6, 10e6, 100e6, 500e6];  % Different sampling rates

figure('Position', [100, 100, 1400, 1000]);

for i = 1:length(fs_values)
    fs = fs_values(i);
    
    % Generate PN code with current fs
    config = gen_pn_code('carrier_freq', carrier_freq, 'fs', fs, ...
                        'modulation', 'sine', 'export_to_base', false);
    
    subplot(2, 2, i);
    
    % Plot first 50 samples to see the pattern
    n_samples = min(50, length(config.time_interp));
    plot(config.time_interp(1:n_samples) * 1e6, ...  % Convert to microseconds
         config.pn_interp_modulated(1:n_samples), 'b-', 'LineWidth', 1.5);
    
    title(sprintf('fs = %.0f MHz', fs/1e6));
    xlabel('Time (μs)');
    ylabel('Amplitude');
    grid on;
    
    % Add text showing key parameters
    chip_duration = 1/fs * 1e6;  % in microseconds
    total_duration = config.total_duration * 1e6;  % in microseconds
    
    text(0.05, 0.95, sprintf('Chip Duration: %.3f μs\nTotal Duration: %.1f μs', ...
         chip_duration, total_duration), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'FontSize', 8);
end

sgtitle('Effect of Sampling Rate on SSTDR Signal Time Scale', 'FontSize', 14);

%% Show the key insight
fprintf('=== SAMPLING RATE EFFECT ANALYSIS ===\n');
fprintf('Carrier Frequency: %.1f MHz (fixed)\n\n', carrier_freq/1e6);
fprintf('%-12s %-15s %-15s %-15s\n', 'Sampling', 'Chip Duration', 'Total Duration', 'Chip Rate');
fprintf('%-12s %-15s %-15s %-15s\n', '(MHz)', '(μs)', '(μs)', '(MHz)');
fprintf('%s\n', repmat('-', 1, 65));

for i = 1:length(fs_values)
    fs = fs_values(i);
    chip_duration = 1/fs * 1e6;  % microseconds
    total_duration = (1023/fs) * 1e6;  % microseconds (1023 chips)
    chip_rate = fs / 1e6;  % MHz
    
    fprintf('%-12.0f %-15.3f %-15.1f %-15.1f\n', ...
            fs/1e6, chip_duration, total_duration, chip_rate);
end

fprintf('\n🔍 KEY INSIGHT:\n');
fprintf('The sampling rate directly controls chip duration!\n');
fprintf('Higher fs → shorter chips → compressed time scale → different graphs\n');
fprintf('\nThis is why your correlation plots change so dramatically!\n'); 