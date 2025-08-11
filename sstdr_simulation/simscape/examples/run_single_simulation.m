%% SIMPLE SSTDR RUN - Minimal Setup
% Just change the MODEL_NAME and hit run!

clear; clc; close all;

%% Setup paths for organized folder structure
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(current_dir);
addpath(fullfile(parent_dir, 'functions', 'simulation'));
addpath(fullfile(parent_dir, 'config'));

%% 🎯 CHANGE THIS TO YOUR MODEL NAME
MODEL_NAME = 'models/sstdr_basic';  % Set to your Simscape model name, or leave empty to use existing data

%% Quick configuration selection - change the number to try different setups
CONFIG_CHOICE = 6;  % Change this number (1-8) to try different configurations

switch CONFIG_CHOICE
    case 1  % Default: 100 kHz chip rate, 100 kHz carrier, 400 kHz sampling
        config = create_simulation_config('chip_rate', 100e3, 'carrier_freq', 100e3, 'fs', 400e3);                    
    case 2  % Fast: 250 kHz chip rate, 250 kHz carrier, 1 MHz sampling
        config = create_simulation_config('chip_rate', 250e3, 'carrier_freq', 250e3, 'fs', 1e6);                  
    case 3  % High Resolution: 200 kHz chip rate, 200 kHz carrier, 800 kHz sampling
        config = create_simulation_config('chip_rate', 200e3, 'carrier_freq', 200e3, 'fs', 800e3);                   
    case 4  % Unmodulated: 125 kHz chip rate, no carrier, 500 kHz sampling
        config = create_simulation_config('chip_rate', 125e3, 'modulation', 'none', 'fs', 500e3);                
    case 5  % Custom Medium: 200 kHz chip rate, 200 kHz carrier, 800 kHz sampling
        config = create_simulation_config('chip_rate', 200e3, 'carrier_freq', 200e3, 'fs', 3200e3);   
    case 6  % Custom High Frequency: 500 kHz chip rate, 500 kHz carrier, 2 MHz sampling (20x longer simulation)
        config = create_simulation_config('chip_rate', 4e6, 'carrier_freq', 4e6, 'fs', 16e6, 'duration', 2e-5, 'positive_only', true);    
    case 7  % Custom Long Range: 50 kHz chip rate, 50 kHz carrier, 200 kHz sampling
        config = create_simulation_config('chip_rate', 50e3, 'carrier_freq', 50e3, 'fs', 200e3);    
    case 8  % Ultra High Resolution: 300 kHz chip rate, unmodulated, 1.2 MHz sampling
        config = create_simulation_config('chip_rate', 300e3, 'modulation', 'none', 'fs', 1.2e6);   
end

%% Run analysis
run_simulation(config, MODEL_NAME, 'verbose', true);

fprintf('\n🎉 Done! Check the plots!\n'); 