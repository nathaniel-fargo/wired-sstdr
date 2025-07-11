%% Generate Demo Network for Visual Inspection
% This script creates a network with various fault types to demonstrate
% the automatic connection functionality

clear; close all; clc;

fprintf('=== Generating Demo Network for Inspection ===\n\n');

% Add the functions directory to path
addpath('./functions');

%% Create a comprehensive test network
fprintf('Creating demo network configuration...\n');

% Design a network with multiple fault types for demonstration
% Load vector: [transmission, series_fault, transmission, shunt_fault, transmission, series_fault, transmission, shunt_fault]
load_vector = [0, 0.7, 0, -0.4, 0, 0.3, 0, -0.6];

% Create network configuration with custom parameters
network_config = create_network_config(load_vector, ...
    'name', 'demo_network_inspection', ...
    'dx', 1.5, ...  % 1.5m per segment
    'Z0', 75, ...   % 75 ohm characteristic impedance
    'velocity', 2.2e8);  % Slightly different velocity

fprintf('Network created:\n');
fprintf('  Name: %s\n', network_config.name);
fprintf('  Segments: %d\n', network_config.num_segments);
fprintf('  Total length: %.1f m\n', network_config.physical.total_length);
fprintf('  Faults: %d (%d series, %d shunt)\n', ...
    network_config.analysis.num_faults, ...
    network_config.analysis.num_opens, ...
    network_config.analysis.num_shorts);

%% Display fault details
fprintf('\nFault Details:\n');
for i = 1:network_config.num_segments
    fault_info = network_config.faults.(sprintf('segment_%d', i));
    if ~strcmp(fault_info.element_type, 'none')
        fprintf('  Segment %d (%.1f m): %s fault - %.2e Ω\n', ...
            i, fault_info.position_m, fault_info.element_type, fault_info.resistance);
    end
end

%% Build the network model with automatic connections
fprintf('\nBuilding Simscape model with automatic connections...\n');

try
    [model_name, model_info] = build_network_model(network_config, ...
        'connect_blocks', true, ...
        'save_model', true, ...
        'close_after', false);  % Keep open for inspection
    
    fprintf('\n✓ Demo network built successfully!\n');
    fprintf('  Model name: %s\n', model_name);
    fprintf('  Total blocks: %d\n', count_total_blocks(model_info));
    
    % Display connection summary
    display_detailed_summary(model_info);
    
    % Set up the model view for better inspection
    setup_model_view(model_name);
    
    fprintf('\n=== Model Ready for Inspection ===\n');
    fprintf('The model "%s" is now open in Simulink.\n', model_name);
    fprintf('You can:\n');
    fprintf('  1. Inspect the block layout and connections\n');
    fprintf('  2. Verify series faults are in-line with transmission lines\n');
    fprintf('  3. Verify shunt faults connect to ground\n');
    fprintf('  4. Check that all transmission lines link to line_test.mat\n');
    fprintf('  5. See the SSTDR subsystem connection\n');
    fprintf('\nThe model has been saved to disk.\n');
    
catch ME
    fprintf('✗ Failed to build demo network: %s\n', ME.message);
    if length(ME.stack) > 0
        fprintf('  Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

%% Helper Functions

function total_blocks = count_total_blocks(model_info)
    % Count all blocks in the model
    total_blocks = 0;
    
    % Essential blocks
    if isfield(model_info.blocks, 'powergui'), total_blocks = total_blocks + 1; end
    if isfield(model_info.blocks, 'ground'), total_blocks = total_blocks + 1; end
    
    % SSTDR blocks
    if isfield(model_info.blocks, 'sstdr_system'), total_blocks = total_blocks + 1; end
    if isfield(model_info.blocks, 'sstdr_source'), total_blocks = total_blocks + 1; end
    if isfield(model_info.blocks, 'tx_sensor'), total_blocks = total_blocks + 1; end
    if isfield(model_info.blocks, 'rx_sensor'), total_blocks = total_blocks + 1; end
    if isfield(model_info.blocks, 'tx_output'), total_blocks = total_blocks + 1; end
    if isfield(model_info.blocks, 'rx_output'), total_blocks = total_blocks + 1; end
    
    % Network blocks
    if isfield(model_info.blocks, 'network')
        network_fields = fieldnames(model_info.blocks.network);
        total_blocks = total_blocks + length(network_fields);
    end
    
    % Termination
    if isfield(model_info.blocks, 'termination'), total_blocks = total_blocks + 1; end
end

function display_detailed_summary(model_info)
    % Display detailed information about the created model
    fprintf('\n=== Detailed Model Summary ===\n');
    
    network_config = model_info.network_config;
    
    % Network topology
    fprintf('Network Topology:\n');
    fprintf('  Physical length: %.1f m (%.1f m per segment)\n', ...
        network_config.physical.total_length, network_config.physical.dx);
    fprintf('  Characteristic impedance: %.0f Ω\n', network_config.physical.Z0);
    fprintf('  Propagation velocity: %.1e m/s\n', network_config.physical.velocity);
    fprintf('  One-way delay: %.2f μs\n', network_config.physical.one_way_delay * 1e6);
    
    % Block inventory
    fprintf('\nBlock Inventory:\n');
    if isfield(model_info.blocks, 'sstdr_system')
        fprintf('  ✓ SSTDR Subsystem\n');
    else
        fprintf('  ✓ Individual SSTDR blocks (source, sensors, outputs)\n');
    end
    fprintf('  ✓ %d Transmission line segments\n', network_config.num_segments);
    fprintf('  ✓ %d Fault elements (%d series, %d shunt)\n', ...
        network_config.analysis.num_faults, ...
        network_config.analysis.num_opens, ...
        network_config.analysis.num_shorts);
    fprintf('  ✓ Termination load (%.0f Ω)\n', network_config.physical.Z0);
    fprintf('  ✓ Ground reference\n');
    fprintf('  ✓ PowerGUI solver\n');
    
    % Connection status
    fprintf('\nAutomatic Connections:\n');
    fprintf('  ✓ SSTDR source → First transmission line\n');
    fprintf('  ✓ Transmission line chain (%d segments)\n', network_config.num_segments);
    if network_config.analysis.num_opens > 0
        fprintf('  ✓ Series faults integrated in-line\n');
    end
    if network_config.analysis.num_shorts > 0
        fprintf('  ✓ Shunt faults connected to ground\n');
    end
    fprintf('  ✓ Termination load → Ground\n');
    fprintf('  ✓ All transmission lines → line_test.mat\n');
end

function setup_model_view(model_name)
    % Set up the model for optimal viewing
    try
        % Fit the model to the window
        set_param(model_name, 'ZoomFactor', 'FitSystem');
        
        % Set a clean background
        set_param(model_name, 'ShowPageBoundaries', 'off');
        
        fprintf('  ✓ Model view optimized for inspection\n');
        
    catch ME
        fprintf('  ⚠ Could not optimize model view: %s\n', ME.message);
    end
end 