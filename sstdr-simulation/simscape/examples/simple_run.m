%% SIMPLE SSTDR RUN - Minimal Setup
% Just change the MODEL_NAME and hit run!

clear; clc; close all;

%% Setup paths for organized folder structure
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'functions'));
addpath(fullfile(parent_dir, 'config'));

%% 🎯 CHANGE THIS TO YOUR MODEL NAME
MODEL_NAME = 'sstdr_basic';  % Set to your Simscape model name, or leave empty to use existing data

%% Quick configuration selection - change the number to try different setups
CONFIG_CHOICE = 4;  % Change this number (1-8) to try different configurations

switch CONFIG_CHOICE
    case 1  % Default
        sstdr_config('default');                    % 100 kHz carrier, 1 MHz sampling
    case 2  % Fast
        sstdr_config('sine_fast');                  % 250 kHz carrier, 2 MHz sampling  
    case 3  % High Resolution
        sstdr_config('high_res');                   % 100 kHz carrier, 5 MHz sampling
    case 4  % Unmodulated
        sstdr_config('unmodulated');                % No carrier, 1 MHz sampling
    case 5  % Custom Medium
        sstdr_custom_config('carrier_freq', 500e3, 'fs', 5e6);   % 500 kHz, 5 MHz
    case 6  % Custom High Frequency
        sstdr_custom_config('carrier_freq', 1e6, 'fs', 10e6);    % 1 MHz, 10 MHz
    case 7  % Custom Long Range
        sstdr_custom_config('carrier_freq', 50e3, 'fs', 2e6);    % 50 kHz, 2 MHz
    case 8  % Ultra High Resolution
        sstdr_custom_config('modulation', 'none', 'fs', 10e6);   % Unmodulated, 10 MHz
end

%% Run analysis
if ~isempty(MODEL_NAME)
    % Configure and run with model
    configure_model_sampling(MODEL_NAME);
    run_sstdr_analysis('skip', MODEL_NAME);
else
    % Use existing simulation data
    cross_correlate('method', 'both');
end

fprintf('\n🎉 Done! Check the plots!\n'); 