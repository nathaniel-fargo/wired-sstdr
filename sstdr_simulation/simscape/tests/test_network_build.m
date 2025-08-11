%% Test Network Building
% Simple test to build a single network model

clear; close all; clc;

%% Add paths
addpath('functions/network_generation');
addpath('config');

%% Create a simple network configuration
% Create a load vector: [transmission, series_fault, transmission, shunt_fault, transmission]
load_vector = [0, 0.5, 0, -0.3, 0];

% Create network configuration
network_config = create_network_config(load_vector, ...
    'dx', 1.0, ...           % 1 meter per segment
    'Z0', 50, ...            % 50 ohm characteristic impedance
    'velocity', 2e8, ...     % 2e8 m/s propagation velocity
    'name', 'test_network');

% Display network info
fprintf('\n=== Network Configuration ===\n');
fprintf('Name: %s\n', network_config.name);
fprintf('Segments: %d\n', network_config.num_segments);
fprintf('Total Length: %.1f m\n', network_config.physical.total_length);
fprintf('Faults: %d\n', network_config.analysis.num_faults);

%% Build the network model
fprintf('\n=== Building Network Model ===\n');

model_name = 'test_network_model';

try
    % Build the network model
    [built_model_name, model_info] = build_network_model(network_config, ...
        'model_name', model_name, ...
        'save_model', false, ...
        'close_after', false, ...
        'connect_blocks', true);
    
    fprintf('\n✓ Network model built successfully: %s\n', built_model_name);
    
    % Display model info
    fprintf('\nModel contains:\n');
    fprintf('  - %d blocks total\n', length(fieldnames(model_info.blocks)));
    fprintf('  - Network blocks: %d\n', length(fieldnames(model_info.blocks.network)));
    
    % Check if model is open
    if bdIsLoaded(model_name)
        fprintf('  - Model is open and ready for inspection\n');
        fprintf('  - You can now view the model in Simulink\n');
    end
    
catch ME
    fprintf('\n✗ Failed to build network model: %s\n', ME.message);
    fprintf('Error details: %s\n', ME.getReport);
end

%% Summary
fprintf('\n=== Test Complete ===\n');
fprintf('If successful, you should see a Simulink model with:\n');
fprintf('  - SSTDR system connected to transmission line chain\n');
fprintf('  - %d transmission line segments\n', network_config.num_segments);
fprintf('  - Series fault at segment 2 (%.1f ohms)\n', network_config.faults.segment_2.resistance);
fprintf('  - Shunt fault at segment 4 (%.1f ohms)\n', network_config.faults.segment_4.resistance);
fprintf('  - Termination load at the end\n'); 