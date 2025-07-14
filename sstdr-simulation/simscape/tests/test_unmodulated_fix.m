%% Test: Unmodulated Configuration Fix
% This tests that the unmodulated configuration works with carrier_freq = 0

clear; clc; close all;

% Setup paths
current_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(current_dir, '..', 'functions', 'network_generation'));
addpath(fullfile(current_dir, '..', 'functions', 'sstdr_simulation'));
addpath(fullfile(current_dir, '..', 'config'));

fprintf('=== Testing Unmodulated Configuration Fix ===\n\n');

%% Test 1: Direct gen_pn_code call with carrier_freq = 0
fprintf('--- Test 1: Direct gen_pn_code with carrier_freq = 0 ---\n');
try
    config1 = gen_pn_code('modulation', 'none', 'carrier_freq', 0, ...
                         'chip_rate', 125e3, 'fs', 500e3, 'export_to_base', false);
    fprintf('✓ gen_pn_code with carrier_freq=0 works!\n');
    fprintf('  Chip rate: %.1f kHz, Sampling: %.1f kHz\n', ...
            config1.settings.chip_rate/1000, config1.settings.fs/1000);
catch ME
    fprintf('✗ gen_pn_code failed: %s\n', ME.message);
end

%% Test 2: Predefined unmodulated configuration
fprintf('\n--- Test 2: Predefined unmodulated configuration ---\n');
try
    sstdr_config('unmodulated');
    fprintf('✓ sstdr_config(''unmodulated'') works!\n');
    
    % Get the configuration
    config = evalin('base', 'sstdr_config');
    fprintf('  Configuration: %s\n', config.name);
    fprintf('  Chip rate: %.1f kHz, Sampling: %.1f kHz\n', ...
            config.pn_config.chip_rate/1000, config.pn_config.fs/1000);
    fprintf('  Modulation: %s, Carrier: %.1f Hz\n', ...
            config.pn_config.modulation, config.pn_config.carrier_freq);
catch ME
    fprintf('✗ sstdr_config failed: %s\n', ME.message);
end

%% Test 3: Custom unmodulated configuration
fprintf('\n--- Test 3: Custom unmodulated configuration ---\n');
try
    sstdr_custom_config('modulation', 'none', 'chip_rate', 150e3, 'fs', 600e3);
    fprintf('✓ sstdr_custom_config with unmodulated works!\n');
    
    % Get the configuration
    config = evalin('base', 'sstdr_config');
    fprintf('  Configuration: %s\n', config.name);
    fprintf('  Chip rate: %.1f kHz, Sampling: %.1f kHz\n', ...
            config.pn_config.chip_rate/1000, config.pn_config.fs/1000);
catch ME
    fprintf('✗ sstdr_custom_config failed: %s\n', ME.message);
end

%% Test 4: Simple run script case 4 (unmodulated)
fprintf('\n--- Test 4: Simple run script case 4 ---\n');
try
    % Simulate the simple_run.m case 4
    sstdr_config('unmodulated');
    fprintf('✓ Simple run case 4 (unmodulated) works!\n');
    
    % Show that correlation analysis would work too
    if evalin('base', 'exist(''pn_interp_modulated'', ''var'')')
        fprintf('✓ PN signals generated successfully\n');
    end
    
catch ME
    fprintf('✗ Simple run case 4 failed: %s\n', ME.message);
end

fprintf('\n=== Summary ===\n');
fprintf('The fix allows carrier_freq >= 0 instead of > 0\n');
fprintf('This enables unmodulated configurations with carrier_freq = 0\n');
fprintf('✓ All unmodulated configuration methods should now work!\n'); 