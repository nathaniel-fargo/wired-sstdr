function [model_name, model_info] = build_simscape_model(load_vector, metadata, varargin)
%BUILD_SIMSCAPE_MODEL Create Simscape model from load vector
%
% Usage:
%   [mdl, info] = build_simscape_model(loads, meta)
%   [mdl, info] = build_simscape_model(loads, meta, 'model_name', 'test_network')
%
% Input:
%   load_vector - Numeric vector of load values from parse_load_string
%   metadata    - Metadata struct from parse_load_string
%   varargin    - Optional name-value pairs:
%                 'model_name' - Name for the model (default: auto-generated)
%                 'save_model' - Save model to disk (default: true)
%                 'close_after' - Close model after building (default: false)
%
% Output:
%   model_name - Name of the created Simulink model
%   model_info - Struct with model details and block information

% Parse optional arguments
p = inputParser;
addParameter(p, 'model_name', '', @ischar);
addParameter(p, 'save_model', true, @islogical);
addParameter(p, 'close_after', false, @islogical);
parse(p, varargin{:});

% Generate model name if not provided
if isempty(p.Results.model_name)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    model_name = sprintf('sstdr_network_%s', timestamp);
else
    model_name = p.Results.model_name;
end

fprintf('Building Simscape model: %s\n', model_name);
fprintf('  Segments: %d, Length: %.1f m\n', metadata.num_segments, metadata.total_length);

% Create new model
try
    new_system(model_name);
    open_system(model_name);
    
    % Set model properties
    set_param(model_name, 'SolverType', 'Variable-step');
    set_param(model_name, 'Solver', 'ode23t');
    set_param(model_name, 'StopTime', '10e-3');
    
    % Initialize model info structure
    model_info = struct();
    model_info.model_name = model_name;
    model_info.num_segments = metadata.num_segments;
    model_info.blocks = struct();
    model_info.connections = {};
    
    % Add Simscape configuration block (using powerlib for compatibility)
    solver_block = add_block('powerlib/powergui', [model_name '/powergui']);
    set_param(solver_block, 'Position', [50, 50, 150, 100]);
    
    % Add electrical reference (ground)
    ground_block = add_block('spsGroundLib/Ground', [model_name '/Ground']);
    set_param(ground_block, 'Position', [50, 200, 100, 230]);
    
    % Build the transmission line network
    [network_blocks, connections] = build_transmission_line_network(model_name, load_vector, metadata);
    
    % Store block information
    model_info.blocks.solver = solver_block;
    model_info.blocks.ground = ground_block;
    model_info.blocks.network = network_blocks;
    model_info.connections = connections;
    
    % Add SSTDR source and measurement points
    [sstdr_blocks, sstdr_connections] = add_sstdr_interface(model_name, network_blocks);
    model_info.blocks.sstdr = sstdr_blocks;
    model_info.connections = [model_info.connections, sstdr_connections];
    
    % Connect all blocks
    connect_all_blocks(model_name, model_info);
    
    % Save model if requested
    if p.Results.save_model
        save_system(model_name);
        fprintf('Model saved: %s.slx\n', model_name);
    end
    
    % Close model if requested
    if p.Results.close_after
        close_system(model_name);
    end
    
    fprintf('Model build complete: %s\n', model_name);
    
catch ME
    % Clean up on error
    if exist(model_name, 'var') && bdIsLoaded(model_name)
        close_system(model_name, 0);
    end
    rethrow(ME);
end

end

function [network_blocks, connections] = build_transmission_line_network(model_name, load_vector, metadata)
%BUILD_TRANSMISSION_LINE_NETWORK Create the main transmission line network

network_blocks = struct();
connections = {};

% Starting position for network blocks
start_x = 200;
start_y = 200;
block_spacing = 150;

% Create transmission line segments and fault elements
for i = 1:length(load_vector)
    segment_name = sprintf('TL_Segment_%d', i);
    
    % Add transmission line block
    tl_block = add_block('spsDistributedParametersLineFrequencyDependentLib/Distributed Parameters Line (Frequency-Dependent)', ...
        [model_name '/' segment_name]);
    
    % Position the block
    x_pos = start_x + (i-1) * block_spacing;
    set_param(tl_block, 'Position', [x_pos, start_y, x_pos+80, start_y+60]);
    
    % Configure transmission line parameters
    set_param(tl_block, 'Z0', sprintf('%.1f', metadata.Z0));
    set_param(tl_block, 'length', sprintf('%.3f', metadata.dx));
    set_param(tl_block, 'velocity', sprintf('%.1e', metadata.velocity));
    
    % Store block reference
    network_blocks.(segment_name) = tl_block;
    
    % Add fault element if needed
    [R_value, element_type] = map_load_to_resistance(load_vector(i));
    
    if ~strcmp(element_type, 'none')
        fault_name = sprintf('Fault_%d_%s', i, element_type);
        fault_block = add_fault_element(model_name, fault_name, element_type, R_value, x_pos, start_y);
        network_blocks.(fault_name) = fault_block;
        
        % Store connection information
        connections{end+1} = struct('type', element_type, 'segment', i, 'resistance', R_value);
    end
end

end

function fault_block = add_fault_element(model_name, fault_name, element_type, R_value, x_pos, y_pos)
%ADD_FAULT_ELEMENT Add series or shunt fault element

if strcmp(element_type, 'series')
    % Add series resistor
    fault_block = add_block('spsSeriesRLCBranchLib/Series RLC Branch', ...
        [model_name '/' fault_name]);
    set_param(fault_block, 'Position', [x_pos, y_pos+80, x_pos+60, y_pos+120]);
    set_param(fault_block, 'R', sprintf('%.1e', R_value));
    set_param(fault_block, 'L', '0');
    set_param(fault_block, 'C', 'inf');
    
elseif strcmp(element_type, 'shunt')
    % Add shunt resistor
    fault_block = add_block('spsSeriesRLCBranchLib/Series RLC Branch', ...
        [model_name '/' fault_name]);
    set_param(fault_block, 'Position', [x_pos, y_pos-80, x_pos+60, y_pos-40]);
    set_param(fault_block, 'R', sprintf('%.1e', R_value));
    set_param(fault_block, 'L', '0');
    set_param(fault_block, 'C', 'inf');
    
else
    error('Unknown fault element type: %s', element_type);
end

end

function [sstdr_blocks, connections] = add_sstdr_interface(model_name, network_blocks)
%ADD_SSTDR_INTERFACE Add SSTDR source and measurement points

sstdr_blocks = struct();
connections = {};

% Add voltage source for SSTDR signal
source_block = add_block('spsControlledVoltageSourceLib/Controlled Voltage Source', ...
    [model_name '/SSTDR_Source']);
set_param(source_block, 'Position', [50, 300, 100, 350]);
sstdr_blocks.source = source_block;

% Add current sensor for transmitted signal
tx_sensor = add_block('spsCurrentMeasurementLib/Current Measurement', ...
    [model_name '/TX_Current_Sensor']);
set_param(tx_sensor, 'Position', [120, 300, 170, 350]);
sstdr_blocks.tx_sensor = tx_sensor;

% Add voltage sensor for received signal
rx_sensor = add_block('spsVoltageMeasurementLib/Voltage Measurement', ...
    [model_name '/RX_Voltage_Sensor']);
set_param(rx_sensor, 'Position', [120, 380, 170, 430]);
sstdr_blocks.rx_sensor = rx_sensor;

% Add output ports for data logging
tx_out = add_block('simulink/Sinks/To Workspace', [model_name '/TX_Output']);
set_param(tx_out, 'Position', [200, 320, 250, 340]);
set_param(tx_out, 'VariableName', 'tx_current');
sstdr_blocks.tx_output = tx_out;

rx_out = add_block('simulink/Sinks/To Workspace', [model_name '/RX_Output']);
set_param(rx_out, 'Position', [200, 400, 250, 420]);
set_param(rx_out, 'VariableName', 'rx_voltage');
sstdr_blocks.rx_output = rx_out;

% Store connection information
connections{end+1} = struct('from', 'source', 'to', 'tx_sensor', 'type', 'electrical');
connections{end+1} = struct('from', 'tx_sensor', 'to', 'network_input', 'type', 'electrical');
connections{end+1} = struct('from', 'rx_sensor', 'to', 'network_input', 'type', 'electrical');
connections{end+1} = struct('from', 'tx_sensor', 'to', 'tx_output', 'type', 'signal');
connections{end+1} = struct('from', 'rx_sensor', 'to', 'rx_output', 'type', 'signal');

end

function connect_all_blocks(model_name, model_info)
%CONNECT_ALL_BLOCKS Connect all blocks in the model

% This is a placeholder - actual connection logic would be more complex
% For now, we'll implement basic connections
fprintf('Connecting blocks in model: %s\n', model_name);

% TODO: Implement actual block connections based on model_info.connections
% This requires careful handling of electrical vs signal connections

end 