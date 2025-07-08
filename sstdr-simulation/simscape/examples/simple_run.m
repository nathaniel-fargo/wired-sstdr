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
CONFIG_CHOICE = 6;  % Change this number (1-8) to try different configurations

switch CONFIG_CHOICE
    case 1  % Default: 100 kHz chip rate, 100 kHz carrier, 400 kHz sampling
        sstdr_config('default');                    
    case 2  % Fast: 250 kHz chip rate, 250 kHz carrier, 1 MHz sampling
        sstdr_config('sine_fast');                  
    case 3  % High Resolution: 200 kHz chip rate, 200 kHz carrier, 800 kHz sampling
        sstdr_config('high_res');                   
    case 4  % Unmodulated: 125 kHz chip rate, no carrier, 500 kHz sampling
        sstdr_config('unmodulated');                
    case 5  % Custom Medium: 200 kHz chip rate, 200 kHz carrier, 800 kHz sampling
        sstdr_custom_config('chip_rate', 200e3, 'carrier_freq', 200e3, 'fs', 3200e3);   
    case 6  % Custom High Frequency: 500 kHz chip rate, 500 kHz carrier, 2 MHz sampling
        sstdr_custom_config('chip_rate', 4e6, 'carrier_freq', 8e6, 'fs', 64e6);    
    case 7  % Custom Long Range: 50 kHz chip rate, 50 kHz carrier, 200 kHz sampling
        sstdr_custom_config('chip_rate', 50e3, 'carrier_freq', 50e3, 'fs', 200e3);    
    case 8  % Ultra High Resolution: 300 kHz chip rate, unmodulated, 1.2 MHz sampling
        sstdr_custom_config('chip_rate', 300e3, 'modulation', 'none', 'fs', 1.2e6);   
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