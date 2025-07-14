%% Test SSTDR Simulation Function
% Simple test of the run_sstdr_simulation function

clear; close all; clc;

%% Add paths
addpath('functions/sstdr_simulation');
addpath('functions/dataset_generation');
addpath('config');

%% Create basic configuration
fprintf('=== Creating Configuration ===\n');

config = create_sstdr_dataset_config( ...
    'chip_rate', 5e6, ...        % 5 MHz chip rate
    'carrier_freq', 5e6, ...     % 5 MHz carrier
    'fs', 20e6, ...             % 20 MHz sampling
    'pn_bits', 11, ...          % 2047 chip sequence
    'duration', 20e-6);         % 20 μs simulation

fprintf('✓ Config: %.1f MHz chip rate, %.1f MHz sampling, %.1f μs duration\n', ...
    config.pn_config.chip_rate/1e6, config.pn_config.fs/1e6, config.simulation_config.stop_time*1e6);

%% Find a model to test with
fprintf('\n=== Finding Model ===\n');

model_candidates = {'sstdr_basic', 'demo_network_inspection'};
test_model = '';

for i = 1:length(model_candidates)
    candidate = model_candidates{i};
    if exist(candidate, 'file') || exist([candidate '.slx'], 'file')
        test_model = candidate;
        fprintf('✓ Found model: %s\n', test_model);
        break;
    end
end

if isempty(test_model)
    fprintf('✗ No model found. Available: %s\n', strjoin(model_candidates, ', '));
    return;
end

%% Run SSTDR simulation
fprintf('\n=== Running SSTDR Simulation ===\n');

try
    [results, sim_data] = run_sstdr_simulation(config, test_model, ...
        'plot_results', true, ...
        'verbose', true);
    
    if results.success
        fprintf('\n✓ SUCCESS!\n');
        fprintf('  Simulation time: %.2f seconds\n', results.sim_time);
        fprintf('  TX samples: %d\n', length(results.tx_signal.Data));
        fprintf('  RX samples: %d\n', length(results.rx_signal.Data));
        
        if isfield(results.correlation, 'peaks') && ~isempty(results.correlation.peaks)
            fprintf('  Peaks found: %d\n', length(results.correlation.peaks.values));
            fprintf('  Peak times: ');
            for i = 1:length(results.correlation.peaks.times)
                fprintf('%.2f μs ', results.correlation.peaks.times(i)*1e6);
            end
            fprintf('\n');
        else
            fprintf('  No peaks detected\n');
        end
    else
        fprintf('✗ FAILED: %s\n', results.error_message);
    end
    
catch ME
    fprintf('✗ ERROR: %s\n', ME.message);
end

fprintf('\n=== Test Complete ===\n'); 