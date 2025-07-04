%% Quick SSTDR Analysis - Just Hit Run!
% Uncomment the configuration you want to test and run this script
% Your Simscape model name (change this to your actual model name)

clear; clc; close all;

%% Setup paths for organized folder structure
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'functions'));
addpath(fullfile(parent_dir, 'config'));

%% CONFIGURATION - Uncomment ONE of these sections
fprintf('=== Quick SSTDR Analysis Script ===\n');

% 🎯 EDIT THIS: Set your Simscape model name
MODEL_NAME = 'sstdr_basic';  % Change this to your actual model name
% If you don't have a model yet, set to '' to skip simulation
% MODEL_NAME = '';  

%% ========================================
%% CONFIGURATION OPTIONS - Uncomment ONE
%% ========================================

%% 1. DEFAULT CONFIGURATION (100 kHz chip rate, 100 kHz carrier, 400 kHz sampling)
% sstdr_config('default');

%% 2. FAST SINE WAVE (250 kHz chip rate, 250 kHz carrier, 1 MHz sampling)
% sstdr_config('sine_fast');

%% 3. HIGH RESOLUTION (200 kHz chip rate, 200 kHz carrier, 800 kHz sampling) 
% sstdr_config('high_res');

%% 4. UNMODULATED (125 kHz chip rate, no carrier, 500 kHz sampling)
% sstdr_config('unmodulated');

%% 5. CUSTOM: 200 kHz chip rate, 200 kHz carrier, 800 kHz sampling
% sstdr_custom_config('chip_rate', 200e3, 'carrier_freq', 200e3, 'fs', 800e3);

%% 6. CUSTOM: 500 kHz chip rate, 500 kHz carrier, 2 MHz sampling (high frequency)
sstdr_custom_config('chip_rate', 500e3, 'carrier_freq', 500e3, 'fs', 2e6);

%% 7. CUSTOM: 50 kHz chip rate, 50 kHz carrier, 200 kHz sampling (long range)
% sstdr_custom_config('chip_rate', 50e3, 'carrier_freq', 50e3, 'fs', 200e3);

%% 8. CUSTOM: 300 kHz chip rate, unmodulated, 1.2 MHz sampling (ultra high resolution)
% sstdr_custom_config('chip_rate', 300e3, 'modulation', 'none', 'fs', 1.2e6);

%% 9. INTERACTIVE MODE - Will prompt you for parameters
% sstdr_custom_config();

%% ========================================
%% AUTOMATIC EXECUTION
%% ========================================

%% Configure model if specified
if ~isempty(MODEL_NAME)
    fprintf('\n--- Configuring Simscape Model ---\n');
    try
        configure_model_sampling(MODEL_NAME);
        fprintf('✓ Model configured successfully\n');
    catch ME
        fprintf('⚠ Model configuration failed: %s\n', ME.message);
        fprintf('  Continuing with manual simulation...\n');
        MODEL_NAME = '';  % Skip automatic simulation
    end
else
    fprintf('\n--- No model specified - Using existing simulation data ---\n');
end

%% Run complete analysis
fprintf('\n--- Running SSTDR Analysis ---\n');
try
    if ~isempty(MODEL_NAME)
        % Run with automatic simulation
        run_sstdr_analysis('skip', MODEL_NAME, 'method', 'both');
    else
        % Use existing simulation data
        cross_correlate('method', 'both');
    end
    
    fprintf('\n✓ SSTDR Analysis Complete!\n');
    
catch ME
    fprintf('\n✗ Analysis failed: %s\n', ME.message);
    fprintf('Make sure you have simulation data in workspace variable "out"\n');
end

%% ========================================
%% QUICK RESULTS SUMMARY
%% ========================================
fprintf('\n=== Quick Results Summary ===\n');

% Display current configuration
if evalin('base', 'exist(''sstdr_config'', ''var'')')
    config = evalin('base', 'sstdr_config');
    fprintf('Configuration: %s\n', config.name);
    fprintf('Modulation: %s\n', config.pn_config.modulation);
    if ~strcmp(config.pn_config.modulation, 'none')
        fprintf('Carrier: %.1f kHz\n', config.pn_config.carrier_freq/1000);
    end
    fprintf('Sampling: %.1f kHz\n', config.pn_config.fs/1000);
end

% Display correlation results if available
if evalin('base', 'exist(''correlation_results'', ''var'')')
    results = evalin('base', 'correlation_results');
    
    if isfield(results, 'peak_times_freq') && ~isempty(results.peak_times_freq)
        fprintf('Detected peaks at: ');
        fprintf('%.1f ', results.peak_times_freq * 1e6);
        fprintf('µs\n');
    end
    
    if isfield(results, 'speedup')
        if results.speedup > 1
            fprintf('Performance: Frequency domain %.1fx faster\n', results.speedup);
        else
            fprintf('Performance: Time domain %.1fx faster\n', 1/results.speedup);
        end
    end
else
    fprintf('No correlation results found\n');
end

%% ========================================
%% ADDITIONAL OPTIONS (Uncomment to use)
%% ========================================

%% OPTION A: Parameter sweep different carrier frequencies
% fprintf('\n--- Parameter Sweep: Carrier Frequencies ---\n');
% carrier_freqs = [100e3, 250e3, 500e3, 1e6];
% for i = 1:length(carrier_freqs)
%     fc = carrier_freqs(i);
%     fprintf('\nTesting %.0f kHz carrier...\n', fc/1000);
%     sstdr_custom_config('carrier_freq', fc, 'fs', 10*fc, 'apply', true);
%     if ~isempty(MODEL_NAME)
%         configure_model_sampling(MODEL_NAME);
%         run_sstdr_analysis('skip', MODEL_NAME, 'method', 'freq', 'save_results', true);
%     end
% end

%% OPTION B: Compare time vs frequency domain performance
% fprintf('\n--- Performance Comparison ---\n');
% if evalin('base', 'exist(''out'', ''var'')')
%     tic; cross_correlate('method', 'time', 'plot_results', false); time_duration = toc;
%     tic; cross_correlate('method', 'freq', 'plot_results', false); freq_duration = toc;
%     fprintf('Time domain: %.2f ms\n', time_duration*1000);
%     fprintf('Frequency domain: %.2f ms\n', freq_duration*1000);
%     fprintf('Speedup: %.1fx\n', time_duration/freq_duration);
% end

%% OPTION C: Save all results to file
% timestamp = datestr(now, 'yyyymmdd_HHMMSS');
% filename = sprintf('sstdr_results_%s.mat', timestamp);
% if evalin('base', 'exist(''correlation_results'', ''var'')')
%     results = evalin('base', 'correlation_results');
%     config = evalin('base', 'sstdr_config');
%     save(filename, 'results', 'config', 'timestamp');
%     fprintf('Results saved to: %s\n', filename);
% end

fprintf('\n=== Script Complete ===\n');
fprintf('Check the figure windows for correlation plots!\n'); 