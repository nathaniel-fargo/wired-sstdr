function [model_name, model_info] = build_simple_model(load_vector, metadata, varargin)
%BUILD_SIMPLE_MODEL Create a basic Simscape model structure
%
% This is a simplified version that focuses on getting the basic structure
% working before adding complex parameter configuration

% Parse optional arguments
p = inputParser;
addParameter(p, 'model_name', '', @ischar);
addParameter(p, 'save_model', false, @islogical);
addParameter(p, 'close_after', false, @islogical);
parse(p, varargin{:});

% Generate model name if not provided
if isempty(p.Results.model_name)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    model_name = sprintf('simple_network_%s', timestamp);
else
    model_name = p.Results.model_name;
end

fprintf('Building simple Simscape model: %s\n', model_name);
fprintf('  Segments: %d, Length: %.1f m\n', metadata.num_segments, metadata.total_length);

% Create new model
try
    new_system(model_name);
    open_system(model_name);
    
    % Set basic model properties
    set_param(model_name, 'SolverType', 'Variable-step');
    set_param(model_name, 'Solver', 'ode23t');
    set_param(model_name, 'StopTime', '10e-3');
    
    % Initialize model info structure
    model_info = struct();
    model_info.model_name = model_name;
    model_info.num_segments = metadata.num_segments;
    model_info.blocks = struct();
    model_info.load_vector = load_vector;
    model_info.metadata = metadata;
    
    % Add basic required blocks
    add_basic_blocks(model_name, model_info);
    
    % Add transmission line segments
    add_transmission_line_segments(model_name, load_vector, metadata, model_info);
    
    % Save model if requested
    if p.Results.save_model
        save_system(model_name);
        fprintf('Model saved: %s.slx\n', model_name);
    end
    
    % Close model if requested
    if p.Results.close_after
        close_system(model_name);
    end
    
    fprintf('Simple model build complete: %s\n', model_name);
    
catch ME
    % Clean up on error
    if exist(model_name, 'var') && bdIsLoaded(model_name)
        close_system(model_name, 0);
    end
    rethrow(ME);
end

end

function add_basic_blocks(model_name, model_info)
%ADD_BASIC_BLOCKS Add the essential blocks every model needs

% Add powergui (Simscape solver configuration)
try
    powergui_block = add_block('powerlib/powergui', [model_name '/powergui']);
    set_param(powergui_block, 'Position', [50, 50, 150, 100]);
    model_info.blocks.powergui = powergui_block;
    fprintf('  ✓ Added powergui block\n');
catch ME
    fprintf('  ✗ Failed to add powergui: %s\n', ME.message);
end

% Add ground reference
try
    ground_block = add_block('spsGroundLib/Ground', [model_name '/Ground']);
    set_param(ground_block, 'Position', [50, 200, 100, 230]);
    model_info.blocks.ground = ground_block;
    fprintf('  ✓ Added ground block\n');
catch ME
    fprintf('  ✗ Failed to add ground: %s\n', ME.message);
end

% Add voltage source
try
    source_block = add_block('spsControlledVoltageSourceLib/Controlled Voltage Source', ...
        [model_name '/SSTDR_Source']);
    set_param(source_block, 'Position', [50, 300, 100, 350]);
    model_info.blocks.source = source_block;
    fprintf('  ✓ Added voltage source\n');
catch ME
    fprintf('  ✗ Failed to add voltage source: %s\n', ME.message);
end

% Add voltage measurement
try
    vmeas_block = add_block('spsVoltageMeasurementLib/Voltage Measurement', ...
        [model_name '/Voltage_Measurement']);
    set_param(vmeas_block, 'Position', [200, 300, 250, 350]);
    model_info.blocks.vmeas = vmeas_block;
    fprintf('  ✓ Added voltage measurement\n');
catch ME
    fprintf('  ✗ Failed to add voltage measurement: %s\n', ME.message);
end

% Add output to workspace
try
    output_block = add_block('simulink/Sinks/To Workspace', [model_name '/Output']);
    set_param(output_block, 'Position', [300, 320, 350, 340]);
    set_param(output_block, 'VariableName', 'sim_output');
    model_info.blocks.output = output_block;
    fprintf('  ✓ Added output block\n');
catch ME
    fprintf('  ✗ Failed to add output: %s\n', ME.message);
end

end

function add_transmission_line_segments(model_name, load_vector, metadata, model_info)
%ADD_TRANSMISSION_LINE_SEGMENTS Add TL segments and fault elements

fprintf('  Adding %d transmission line segments...\n', length(load_vector));

% Starting position for network blocks
start_x = 150;
start_y = 400;
block_spacing = 120;

% Initialize network blocks structure
if ~isfield(model_info.blocks, 'network')
    model_info.blocks.network = struct();
end

for i = 1:length(load_vector)
    segment_name = sprintf('TL_Segment_%d', i);
    
    try
        % Add transmission line block (don't set parameters yet)
        tl_block = add_block('spsDistributedParametersLineFrequencyDependentLib/Distributed Parameters Line (Frequency-Dependent)', ...
            [model_name '/' segment_name]);
        
        % Position the block
        x_pos = start_x + (i-1) * block_spacing;
        set_param(tl_block, 'Position', [x_pos, start_y, x_pos+80, start_y+60]);
        
        % Only set the length parameter for now
        try
            set_param(tl_block, 'L', sprintf('%.3f', metadata.dx));
            fprintf('    ✓ Segment %d: length set to %.3f m\n', i, metadata.dx);
        catch
            fprintf('    ⚠ Segment %d: could not set length parameter\n', i);
        end
        
        % Store block reference
        model_info.blocks.network.(segment_name) = tl_block;
        
        % Add fault element if needed
        [R_value, element_type] = map_load_to_resistance(load_vector(i));
        
        if ~strcmp(element_type, 'none')
            fault_name = sprintf('Fault_%d_%s', i, element_type);
            try
                fault_block = add_block('spsSeriesRLCBranchLib/Series RLC Branch', ...
                    [model_name '/' fault_name]);
                
                % Position fault element
                if strcmp(element_type, 'series')
                    set_param(fault_block, 'Position', [x_pos, start_y+80, x_pos+60, start_y+120]);
                else % shunt
                    set_param(fault_block, 'Position', [x_pos, start_y-80, x_pos+60, start_y-40]);
                end
                
                % Configure as pure resistor
                set_param(fault_block, 'Resistance', sprintf('%.1e', R_value));
                set_param(fault_block, 'Inductance', '0');
                set_param(fault_block, 'Capacitance', 'inf');
                
                model_info.blocks.network.(fault_name) = fault_block;
                fprintf('    ✓ Fault %d: %s %.1e Ω\n', i, element_type, R_value);
                
            catch ME
                fprintf('    ✗ Failed to add fault %d: %s\n', i, ME.message);
            end
        end
        
    catch ME
        fprintf('    ✗ Failed to add segment %d: %s\n', i, ME.message);
    end
end

fprintf('  Network structure complete\n');

end 