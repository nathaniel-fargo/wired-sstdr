%% Test Script for Simple Model Builder
% This script tests the basic model building functionality

clear; clc; close all;

fprintf('=== Testing Simple Model Builder ===\n\n');

% Add necessary paths
addpath('functions');
addpath('config');

%% Test: Build a simple model
fprintf('Building simple model with basic network...\n');

% Create a simple test network
test_string = '(0)(0.5)(-0.3)(0)';
[loads, meta] = parse_load_string(test_string);

try
    % Build the model
    [model_name, model_info] = build_simple_model(loads, meta, ...
        'model_name', 'test_simple_network', ...
        'save_model', true, ...
        'close_after', false);
    
    fprintf('\n✓ Model built successfully: %s\n', model_name);
    if isfield(model_info.blocks, 'network')
        fprintf('  Model has %d network blocks\n', length(fieldnames(model_info.blocks.network)));
    else
        fprintf('  No network blocks field found\n');
    end
    
    % Show what was created
    fprintf('\nBlocks created:\n');
    if isfield(model_info.blocks, 'powergui')
        fprintf('  ✓ powergui (solver)\n');
    end
    if isfield(model_info.blocks, 'ground')
        fprintf('  ✓ ground (reference)\n');
    end
    if isfield(model_info.blocks, 'source')
        fprintf('  ✓ voltage source\n');
    end
    if isfield(model_info.blocks, 'vmeas')
        fprintf('  ✓ voltage measurement\n');
    end
    if isfield(model_info.blocks, 'output')
        fprintf('  ✓ output to workspace\n');
    end
    
    % Show network structure
    fprintf('\nNetwork structure:\n');
    if isfield(model_info.blocks, 'network')
        network_fields = fieldnames(model_info.blocks.network);
        for i = 1:length(network_fields)
            fprintf('  - %s\n', network_fields{i});
        end
    else
        fprintf('  (Network blocks not accessible from returned structure)\n');
    end
    
    fprintf('\nModel is open in Simulink. You can now:\n');
    fprintf('1. Examine the block layout\n');
    fprintf('2. Check block parameters\n');
    fprintf('3. Manually connect blocks if needed\n');
    fprintf('4. Close with: close_system(''%s'', 0)\n', model_name);
    
catch ME
    fprintf('✗ Model building failed: %s\n', ME.message);
    if length(ME.stack) > 0
        fprintf('  Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

fprintf('\n=== Test Complete ===\n'); 