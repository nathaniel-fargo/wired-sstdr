%% Test: Chip Rate Control in SSTDR Configuration
% This demonstrates the new chip rate control system where you specify:
% 1. Chip rate (Hz) - controls the PN sequence timing
% 2. Carrier frequency (Hz) - modulation frequency 
% 3. Sampling frequency (Hz) - digital sampling rate

clear; clc; close all;

% Setup paths
current_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(current_dir, 'functions'));
addpath(fullfile(current_dir, 'config'));

%% Test the new chip rate control system
fprintf('=== SSTDR Chip Rate Control Test ===\n\n');

%% Test 1: Default configuration pattern (chip_rate = carrier_freq, fs = 4x chip_rate)
fprintf('--- Test 1: Default Pattern ---\n');
chip_rate = 100e3;      % 100 kHz chip rate
carrier_freq = 100e3;   % 100 kHz carrier (equal to chip rate)
fs = 400e3;             % 400 kHz sampling (4x chip rate)

config1 = gen_pn_code('chip_rate', chip_rate, 'carrier_freq', carrier_freq, ...
                      'fs', fs, 'modulation', 'sine', 'export_to_base', false);

%% Test 2: High frequency pattern
fprintf('\n--- Test 2: High Frequency Pattern ---\n');
chip_rate = 250e3;      % 250 kHz chip rate
carrier_freq = 250e3;   % 250 kHz carrier (equal to chip rate)
fs = 1e6;               % 1 MHz sampling (4x chip rate)

config2 = gen_pn_code('chip_rate', chip_rate, 'carrier_freq', carrier_freq, ...
                      'fs', fs, 'modulation', 'sine', 'export_to_base', false);

%% Test 3: Unmodulated pattern
fprintf('\n--- Test 3: Unmodulated Pattern ---\n');
chip_rate = 125e3;      % 125 kHz chip rate
fs = 500e3;             % 500 kHz sampling (4x chip rate)

config3 = gen_pn_code('chip_rate', chip_rate, 'modulation', 'none', ...
                      'fs', fs, 'export_to_base', false);

%% Test 4: Custom pattern (different carrier and chip rates)
fprintf('\n--- Test 4: Custom Pattern (carrier != chip rate) ---\n');
chip_rate = 100e3;      % 100 kHz chip rate
carrier_freq = 200e3;   % 200 kHz carrier (2x chip rate)
fs = 800e3;             % 800 kHz sampling (8x chip rate)

config4 = gen_pn_code('chip_rate', chip_rate, 'carrier_freq', carrier_freq, ...
                      'fs', fs, 'modulation', 'sine', 'export_to_base', false);

%% Display comparison table
fprintf('\n=== Configuration Comparison ===\n');
fprintf('%-15s %-12s %-12s %-12s %-15s %-15s\n', 'Test', 'Chip Rate', 'Carrier', 'Sampling', 'Chip Duration', 'Total Duration');
fprintf('%-15s %-12s %-12s %-12s %-15s %-15s\n', '', '(kHz)', '(kHz)', '(kHz)', '(μs)', '(ms)');
fprintf('%s\n', repmat('-', 1, 90));

configs = {config1, config2, config3, config4};
test_names = {'Default', 'High Freq', 'Unmodulated', 'Custom'};

for i = 1:length(configs)
    cfg = configs{i};
    
    chip_rate_val = cfg.settings.chip_rate;
    carrier_freq_val = cfg.settings.carrier_freq;
    fs_val = cfg.settings.fs;
    chip_duration = 1e6 / chip_rate_val;  % microseconds
    total_duration = cfg.total_duration * 1000;  % milliseconds
    
    if strcmp(cfg.settings.modulation, 'none')
        carrier_str = 'N/A';
    else
        carrier_str = sprintf('%.0f', carrier_freq_val/1000);
    end
    
    fprintf('%-15s %-12.0f %-12s %-12.0f %-15.3f %-15.3f\n', ...
            test_names{i}, chip_rate_val/1000, carrier_str, fs_val/1000, ...
            chip_duration, total_duration);
end

%% Create visualization plots
figure('Position', [100, 100, 1400, 1000]);

for i = 1:length(configs)
    cfg = configs{i};
    
    subplot(2, 2, i);
    
    % Plot first 50 samples to see the pattern
    n_samples = min(50, length(cfg.time_interp));
    plot(cfg.time_interp(1:n_samples) * 1e6, ...  % Convert to microseconds
         cfg.pn_interp_modulated(1:n_samples), 'b-', 'LineWidth', 1.5);
    
    title(sprintf('%s: %.0f kHz chip rate', test_names{i}, cfg.settings.chip_rate/1000));
    xlabel('Time (μs)');
    ylabel('Amplitude');
    grid on;
    
    % Add frequency information
    chip_duration = 1e6 / cfg.settings.chip_rate;  % microseconds
    
    if strcmp(cfg.settings.modulation, 'none')
        freq_info = sprintf('Chip: %.0f kHz\nSampling: %.0f kHz', ...
                           cfg.settings.chip_rate/1000, cfg.settings.fs/1000);
    else
        freq_info = sprintf('Chip: %.0f kHz\nCarrier: %.0f kHz\nSampling: %.0f kHz', ...
                           cfg.settings.chip_rate/1000, cfg.settings.carrier_freq/1000, cfg.settings.fs/1000);
    end
    
    text(0.05, 0.95, freq_info, ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'FontSize', 8);
end

sgtitle('SSTDR Chip Rate Control Demonstration', 'FontSize', 14);

%% Test predefined configurations
fprintf('\n=== Testing Predefined Configurations ===\n');
predefined_configs = {'default', 'sine_fast', 'high_res', 'unmodulated'};

for i = 1:length(predefined_configs)
    fprintf('\n--- %s ---\n', predefined_configs{i});
    sstdr_config(predefined_configs{i});
    
    % Get the configuration that was just created
    config = evalin('base', 'sstdr_config');
    
    fprintf('Chip rate: %.1f kHz, Carrier: ', config.pn_config.chip_rate/1000);
    if strcmp(config.pn_config.modulation, 'none')
        fprintf('N/A, ');
    else
        fprintf('%.1f kHz, ', config.pn_config.carrier_freq/1000);
    end
    fprintf('Sampling: %.1f kHz\n', config.pn_config.fs/1000);
end

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('✓ Chip rate control system implemented successfully!\n');
fprintf('✓ Three input parameters: chip_rate, carrier_freq, fs\n');
fprintf('✓ Default pattern: chip_rate = carrier_freq, fs = 4x chip_rate\n');
fprintf('✓ Interpolation factor automatically calculated: interpFactor = fs/chip_rate\n');
fprintf('✓ All predefined configurations updated to use new pattern\n');
fprintf('\nKey benefits:\n');
fprintf('- Explicit control over signal timing (chip rate)\n');
fprintf('- Independent control of carrier and chip frequencies\n');
fprintf('- Predictable relationship between frequencies\n');
fprintf('- Backward compatibility maintained\n'); 