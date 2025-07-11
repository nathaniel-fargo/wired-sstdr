%% Test Script for Network Builder Improvements
% This script tests the improvements: series fault placement, SSTDR subsystem, and line_test.mat linking

clear; clc; close all;

fprintf('=== Testing Network Builder Improvements ===\n\n');

% Add necessary paths
addpath('functions');
addpath('config');

%% Test 1: Create a network with both series and shunt faults
fprintf('Test 1: Creating network with mixed fault types\n');

try
    % Create a network with specific fault pattern
    load_vector = [0, 0.6, -0.4, 0, 0.8, -0.2];  % Series, shunt, series, shunt
    test_config = create_network_config(load_vector, 'name', 'mixed_faults_test');
    
    fprintf('  ✓ Created network with %d faults\n', test_config.analysis.num_faults);
    fprintf('    - Series faults: %d\n', test_config.analysis.num_opens);
    fprintf('    - Shunt faults: %d\n', test_config.analysis.num_shorts);
    
    % Show fault details
    fprintf('  Fault details:\n');
    for i = 1:test_config.num_segments
        fault_info = test_config.faults.(sprintf('segment_%d', i));
        if ~strcmp(fault_info.element_type, 'none')
            fprintf('    Segment %d: %s fault, %.1e Ω\n', ...
                i, fault_info.element_type, fault_info.resistance);
        end
    end
    
catch ME
    fprintf('  ✗ Failed to create network: %s\n', ME.message);
end

%% Test 2: Build model with improved block placement
fprintf('\nTest 2: Building model with improved block placement\n');

try
    if exist('test_config', 'var')
        [model_name, model_info] = build_network_model(test_config, ...
            'model_name', 'test_improvements', ...
            'save_model', true, ...
            'close_after', false);
        
        fprintf('  ✓ Model built successfully: %s\n', model_name);
        
        % Check what SSTDR blocks were created
        if isfield(model_info.blocks, 'sstdr_system')
            fprintf('  ✓ SSTDR subsystem block added\n');
        elseif isfield(model_info.blocks, 'sstdr_source')
            fprintf('  ⚠ Fell back to individual SSTDR blocks\n');
        else
            fprintf('  ✗ No SSTDR blocks found\n');
        end
        
        % Check network structure
        if isfield(model_info.blocks, 'network')
            network_blocks = fieldnames(model_info.blocks.network);
            tl_blocks = network_blocks(contains(network_blocks, 'TL_Segment'));
            fault_blocks = network_blocks(contains(network_blocks, 'Fault'));
            
            fprintf('  Network structure:\n');
            fprintf('    - Transmission line segments: %d\n', length(tl_blocks));
            fprintf('    - Fault elements: %d\n', length(fault_blocks));
            
            % Show fault block types
            for i = 1:length(fault_blocks)
                if contains(fault_blocks{i}, 'series')
                    fprintf('      - %s (in-line)\n', fault_blocks{i});
                else
                    fprintf('      - %s (to ground)\n', fault_blocks{i});
                end
            end
        end
        
    else
        fprintf('  ✗ No test config available\n');
    end
    
catch ME
    fprintf('  ✗ Model building failed: %s\n', ME.message);
    if length(ME.stack) > 0
        fprintf('    Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

%% Test 3: Check line_test.mat file
fprintf('\nTest 3: Checking line_test.mat file\n');

try
    if exist('line_test.mat', 'file')
        fprintf('  ✓ line_test.mat file found\n');
        
        % Load and inspect the file
        line_data = load('line_test.mat');
        fprintf('  File contents:\n');
        fields = fieldnames(line_data);
        for i = 1:length(fields)
            fprintf('    - %s: %s\n', fields{i}, class(line_data.(fields{i})));
        end
        
    else
        fprintf('  ✗ line_test.mat file not found\n');
        fprintf('  ⚠ Transmission lines may not have proper frequency data\n');
    end
    
catch ME
    fprintf('  ✗ Error checking line_test.mat: %s\n', ME.message);
end

%% Summary and recommendations
fprintf('\n=== Test Summary ===\n');

if exist('model_name', 'var') && bdIsLoaded(model_name)
    fprintf('✓ Model "%s" is open in Simulink\n', model_name);
    fprintf('  Recommendations:\n');
    fprintf('  1. Examine series fault placement - should be in-line with TL segments\n');
    fprintf('  2. Check shunt fault placement - should connect to ground\n');
    fprintf('  3. Verify transmission line parameters link to line_test.mat\n');
    fprintf('  4. Create a custom SSTDR subsystem block for reuse\n');
    fprintf('  5. Test manual connections between blocks\n');
    
    fprintf('\n  Close model with: close_system(''%s'', 0)\n', model_name);
else
    fprintf('✗ No model was created successfully\n');
end

fprintf('\n=== Next Steps ===\n');
fprintf('1. Create SSTDR subsystem block in sstdr_basic model\n');
fprintf('2. Implement automatic block connections\n');
fprintf('3. Test end-to-end SSTDR simulation\n');
fprintf('4. Verify frequency-dependent transmission line behavior\n');

fprintf('\n=== Test Complete ===\n'); 