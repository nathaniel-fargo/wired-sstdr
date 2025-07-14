%% Test Script for Network Builder Functions
% This script tests the basic functionality of the network builder

clear; clc; close all;

fprintf('=== Testing Network Builder Functions ===\n\n');

%% Test 1: Parse Load String
fprintf('Test 1: Parse Load String\n');
test_string = '(0)(0.3)(-0.5)(1)(0)(-1)';
try
    [loads, meta] = parse_load_string(test_string);
    fprintf('✓ Successfully parsed: %s\n', test_string);
    fprintf('  Load vector: [%s]\n', sprintf('%.1f ', loads));
    fprintf('  Total length: %.1f m\n', meta.total_length);
    fprintf('  One-way delay: %.1f μs\n', meta.one_way_delay * 1e6);
catch ME
    fprintf('✗ Parse failed: %s\n', ME.message);
end

%% Test 2: Resistance Mapping
fprintf('\nTest 2: Resistance Mapping\n');
test_values = [-1, -0.5, 0, 0.5, 1];
for val = test_values
    try
        [R, type] = map_load_to_resistance(val);
        if strcmp(type, 'none')
            fprintf('  Load %.1f → %s (no element)\n', val, type);
        else
            fprintf('  Load %.1f → %s %.1e Ω\n', val, type, R);
        end
    catch ME
        fprintf('  ✗ Mapping failed for %.1f: %s\n', val, ME.message);
    end
end

%% Test 3: Model Building (Basic Structure)
fprintf('\nTest 3: Model Building\n');
try
    % Use a simple test case
    simple_string = '(0)(0.3)(-0.5)';
    [loads, meta] = parse_load_string(simple_string);
    
    % Build model (but don't save or keep open)
    [model_name, model_info] = build_simscape_model(loads, meta, ...
        'model_name', 'test_network_basic', ...
        'save_model', false, ...
        'close_after', true);
    
    fprintf('✓ Model built successfully: %s\n', model_name);
    fprintf('  Segments: %d\n', model_info.num_segments);
    fprintf('  Network blocks: %d\n', length(fieldnames(model_info.blocks.network)));
    
catch ME
    fprintf('✗ Model building failed: %s\n', ME.message);
    fprintf('  Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

%% Test 4: Integration with SSTDR Config
fprintf('\nTest 4: Integration with SSTDR Config\n');
try
    % Load SSTDR configuration
    config = sstdr_config('default');
    fprintf('✓ SSTDR config loaded: %s\n', config.name);
    fprintf('  Chip rate: %.1f kHz\n', config.pn_config.chip_rate/1000);
    fprintf('  Sampling rate: %.1f kHz\n', config.pn_config.fs/1000);
    
    % Show how the two would integrate
    fprintf('  Simulation time: %.1f ms\n', config.simulation_config.stop_time*1000);
    fprintf('  Max detectable distance: %.1f m\n', ...
        config.simulation_config.stop_time * meta.velocity / 2);
    
catch ME
    fprintf('✗ SSTDR config failed: %s\n', ME.message);
end

fprintf('\n=== Test Complete ===\n');

%% Helper: Show what still needs to be implemented
fprintf('\n=== Next Steps for Full Implementation ===\n');
fprintf('1. Complete the connect_all_blocks() function\n');
fprintf('2. Add proper SSTDR signal generation (from existing gen_pn_code)\n');
fprintf('3. Connect to run_sstdr_analysis.m for end-to-end testing\n');
fprintf('4. Add termination resistor at far end of network\n');
fprintf('5. Implement proper electrical connections between TL segments\n');
fprintf('6. Add error handling for invalid Simscape block paths\n');
fprintf('7. Test with actual Simscape simulation\n'); 