function [model_name, model_info] = build_network_model(network_config, varargin)
%BUILD_NETWORK_MODEL Create Simscape model from network configuration structure
%
% Usage:
%   [mdl, info] = build_network_model(network_config)
%   [mdl, info] = build_network_model(network_config, 'save_model', true)
%
% Input:
%   network_config - Network configuration structure from create_network_config
%   varargin       - Optional name-value pairs:
%                    'model_name' - Name for the model (default: from config)
%                    'save_model' - Save model to disk (default: false)
%                    'close_after' - Close model after building (default: false)
%                    'connect_blocks' - Automatically connect blocks (default: false)
%
% Output:
%   model_name - Name of the created Simulink model
%   model_info - Struct with model details and block information

% Parse optional arguments
p = inputParser;
addRequired(p, 'network_config', @isstruct);
addParameter(p, 'model_name', '', @ischar);
addParameter(p, 'save_model', false, @islogical);
addParameter(p, 'close_after', false, @islogical);
addParameter(p, 'connect_blocks', false, @islogical);
parse(p, network_config, varargin{:});

% Use network name or generate model name
if isempty(p.Results.model_name)
    model_name = network_config.name;
else
    model_name = p.Results.model_name;
end

fprintf('Building Simscape model: %s\n', model_name);
fprintf('  Network: %s (%d segments)\n', network_config.name, network_config.num_segments);
fprintf('  Length: %.1f m, Faults: %d\n', network_config.physical.total_length, network_config.analysis.num_faults);

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
    model_info.network_config = network_config;
    model_info.blocks = struct();
    model_info.connections = {};
    
    % Add essential blocks
    model_info = add_essential_blocks(model_name, model_info);
    
    % Add SSTDR interface blocks
    model_info = add_sstdr_interface_blocks(model_name, model_info);
    
    % Add transmission line network
    model_info = add_transmission_line_network(model_name, network_config, model_info);
    
    % Connect blocks if requested
    if p.Results.connect_blocks
        connect_network_blocks(model_name, model_info);
    end
    
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
    if exist('model_name', 'var') && bdIsLoaded(model_name)
        close_system(model_name, 0);
    end
    rethrow(ME);
end

end

function model_info = add_essential_blocks(model_name, model_info)
%ADD_ESSENTIAL_BLOCKS Add blocks required for any Simscape model

% Add powergui (Simscape solver configuration)
try
    powergui_block = add_block('powerlib/powergui', [model_name '/powergui']);
    set_param(powergui_block, 'Position', [50, 50, 150, 100]);
    
    % Set to discrete mode instead of continuous
    set_param(powergui_block, 'SimulationMode', 'Discrete');
    set_param(powergui_block, 'SampleTime', 'sim_max_step');  % Use workspace variable
    
    model_info.blocks.powergui = powergui_block;
    fprintf('  ✓ Added powergui block (discrete mode, sample time: sim_max_step)\n');
catch ME
    fprintf('  ✗ Failed to add powergui: %s\n', ME.message);
end

% Add main ground reference (will be positioned at end of line later)
try
    ground_block = add_block('spsGroundLib/Ground', [model_name '/Ground']);
    % Position will be set after we know the network length
    model_info.blocks.ground = ground_block;
    fprintf('  ✓ Added ground block\n');
catch ME
    fprintf('  ✗ Failed to add ground: %s\n', ME.message);
end

end

function model_info = add_sstdr_interface_blocks(model_name, model_info)
%ADD_SSTDR_INTERFACE_BLOCKS Add SSTDR subsystem block

% Try to add custom SSTDR subsystem block first
try
    % First try to use a custom SSTDR subsystem if it exists
    sstdr_block = add_block('sstdr_basic/SSTDR_Subsystem', [model_name '/SSTDR_System']);
    set_param(sstdr_block, 'Position', [50, 250, 150, 350]);
    model_info.blocks.sstdr_system = sstdr_block;
    fprintf('  ✓ Added SSTDR subsystem block\n');
    
    % Add dedicated ground for SSTDR subsystem
    sstdr_ground = add_block('spsGroundLib/Ground', [model_name '/SSTDR_Ground']);
    set_param(sstdr_ground, 'Position', [50, 370, 100, 400]);
    model_info.blocks.sstdr_ground = sstdr_ground;
    fprintf('  ✓ Added SSTDR ground block\n');
    
catch ME1
    % If custom subsystem doesn't exist, try to create a subsystem from the existing model
    try
        % Create a subsystem that references the existing SSTDR blocks
        sstdr_subsys = add_block('simulink/Ports & Subsystems/Subsystem', [model_name '/SSTDR_System']);
        set_param(sstdr_subsys, 'Position', [50, 250, 150, 350]);
        model_info.blocks.sstdr_system = sstdr_subsys;
        fprintf('  ✓ Added SSTDR subsystem placeholder\n');
        fprintf('  ⚠ Custom SSTDR block not found, using placeholder\n');
        fprintf('  ⚠ You should create an SSTDR subsystem block for reuse\n');
        
        % Add dedicated ground for SSTDR subsystem
        sstdr_ground = add_block('spsGroundLib/Ground', [model_name '/SSTDR_Ground']);
        set_param(sstdr_ground, 'Position', [50, 370, 100, 400]);
        model_info.blocks.sstdr_ground = sstdr_ground;
        fprintf('  ✓ Added SSTDR ground block\n');
        
    catch ME2
        % Fall back to individual blocks
        fprintf('  ⚠ Could not add SSTDR subsystem, falling back to individual blocks\n');
        model_info = add_individual_sstdr_blocks(model_name, model_info);
    end
end

end

function model_info = add_individual_sstdr_blocks(model_name, model_info)
%ADD_INDIVIDUAL_SSTDR_BLOCKS Fallback to individual SSTDR blocks

% Add controlled voltage source for SSTDR signal
try
    source_block = add_block('spsControlledVoltageSourceLib/Controlled Voltage Source', ...
        [model_name '/SSTDR_Source']);
    set_param(source_block, 'Position', [50, 250, 100, 300]);
    model_info.blocks.sstdr_source = source_block;
    fprintf('  ✓ Added SSTDR voltage source\n');
catch ME
    fprintf('  ✗ Failed to add SSTDR source: %s\n', ME.message);
end

% Add current measurement for transmitted signal
try
    tx_sensor = add_block('spsCurrentMeasurementLib/Current Measurement', ...
        [model_name '/TX_Current_Sensor']);
    set_param(tx_sensor, 'Position', [150, 250, 200, 300]);
    model_info.blocks.tx_sensor = tx_sensor;
    fprintf('  ✓ Added TX current sensor\n');
catch ME
    fprintf('  ✗ Failed to add TX sensor: %s\n', ME.message);
end

% Add voltage measurement for received signal
try
    rx_sensor = add_block('spsVoltageMeasurementLib/Voltage Measurement', ...
        [model_name '/RX_Voltage_Sensor']);
    set_param(rx_sensor, 'Position', [150, 320, 200, 370]);
    model_info.blocks.rx_sensor = rx_sensor;
    fprintf('  ✓ Added RX voltage sensor\n');
catch ME
    fprintf('  ✗ Failed to add RX sensor: %s\n', ME.message);
end

% Add output blocks
try
    tx_out = add_block('simulink/Sinks/To Workspace', [model_name '/TX_Output']);
    set_param(tx_out, 'Position', [250, 260, 300, 280]);
    set_param(tx_out, 'VariableName', 'tx_current');
    model_info.blocks.tx_output = tx_out;
    
    rx_out = add_block('simulink/Sinks/To Workspace', [model_name '/RX_Output']);
    set_param(rx_out, 'Position', [250, 340, 300, 360]);
    set_param(rx_out, 'VariableName', 'rx_voltage');
    model_info.blocks.rx_output = rx_out;
    
    fprintf('  ✓ Added output blocks\n');
catch ME
    fprintf('  ✗ Failed to add output blocks: %s\n', ME.message);
end

end

function model_info = add_transmission_line_network(model_name, network_config, model_info)
%ADD_TRANSMISSION_LINE_NETWORK Add transmission line segments and fault elements

fprintf('  Adding transmission line network...\n');

% Starting position for network blocks
start_x = 350;
start_y = 275;
block_spacing = 180;  % Increased spacing to fit series faults better

% Initialize network blocks structure
model_info.blocks.network = struct();

% Add each segment
for i = 1:network_config.num_segments
    segment_name = sprintf('TL_Segment_%d', i);
    
    try
        % Add transmission line block
        tl_block = add_block('spsDistributedParametersLineFrequencyDependentLib/Distributed Parameters Line (Frequency-Dependent)', ...
            [model_name '/' segment_name]);
        
        % Position the block
        x_pos = start_x + (i-1) * block_spacing;
        set_param(tl_block, 'Position', [x_pos, start_y, x_pos+80, start_y+60]);
        
        % Set transmission line parameters
        set_param(tl_block, 'L', sprintf('%.3f', network_config.physical.dx));
        
        % Link to transmission line file for frequency-dependent parameters
        tl_file = 'line_test.mat';  % Default file
        if isfield(network_config, 'transmission_line_file')
            tl_file = network_config.transmission_line_file;
        elseif isfield(network_config, 'network') && isfield(network_config.network, 'transmission_line_file')
            tl_file = network_config.network.transmission_line_file;
        end
        
        % Try different paths for transmission line file
        tl_paths = {
            fullfile('config', 'line_params', tl_file),  % New location
            fullfile('config', tl_file),                 % Config directory
            tl_file                                      % Current directory
        };
        
        linked = false;
        for path_idx = 1:length(tl_paths)
            try
                set_param(tl_block, 'WBfile', tl_paths{path_idx});
                fprintf('      ✓ Linked to %s\n', tl_paths{path_idx});
                linked = true;
                break;
            catch ME_mat
                % Continue to next path
            end
        end
        
        if ~linked
            fprintf('      ⚠ Could not link to %s in any location\n', tl_file);
        end
        
        % Store block reference
        model_info.blocks.network.(segment_name) = tl_block;
        
        fprintf('    ✓ Segment %d: %.3f m\n', i, network_config.physical.dx);
        
        % Add fault element if needed
        fault_info = network_config.faults.(sprintf('segment_%d', i));
        
        if ~strcmp(fault_info.element_type, 'none')
            fault_name = sprintf('Fault_%d_%s', i, fault_info.element_type);
            
            try
                fault_block = add_block('spsSeriesRLCBranchLib/Series RLC Branch', ...
                    [model_name '/' fault_name]);
                
                % Position fault element
                if strcmp(fault_info.element_type, 'series')
                    % Series fault goes in-line with transmission line (centered in spacing)
                    set_param(fault_block, 'Position', [x_pos+90, start_y+10, x_pos+150, start_y+50]);
                else % shunt
                    % Shunt fault goes to ground (below transmission line)
                    set_param(fault_block, 'Position', [x_pos+20, start_y+80, x_pos+80, start_y+120]);
                    
                    % Add individual ground for this shunt fault
                    ground_name = sprintf('Ground_%d', i);
                    try
                        shunt_ground = add_block('spsGroundLib/Ground', [model_name '/' ground_name]);
                        set_param(shunt_ground, 'Position', [x_pos+20, start_y+140, x_pos+70, start_y+170]);
                        model_info.blocks.network.(ground_name) = shunt_ground;
                        fprintf('      ✓ Added individual ground for shunt fault\n');
                    catch ME_ground
                        fprintf('      ⚠ Could not add individual ground: %s\n', ME_ground.message);
                    end
                end
                
                % Configure as pure resistor (R-only configuration)
                set_param(fault_block, 'Resistance', sprintf('%.1e', fault_info.resistance));
                set_param(fault_block, 'Inductance', '0');
                set_param(fault_block, 'Capacitance', 'inf');
                set_param(fault_block, 'BranchType', 'R');  % Set to display as R instead of RLC
                
                model_info.blocks.network.(fault_name) = fault_block;
                fprintf('      ✓ Fault: %s %.1e Ω\n', fault_info.element_type, fault_info.resistance);
                
            catch ME
                fprintf('      ✗ Failed to add fault: %s\n', ME.message);
            end
        end
        
    catch ME
        fprintf('    ✗ Failed to add segment %d: %s\n', i, ME.message);
    end
end

% Add termination load
try
    term_block = add_block('spsSeriesRLCBranchLib/Series RLC Branch', ...
        [model_name '/Termination_Load']);
    
    x_pos = start_x + network_config.num_segments * block_spacing;
    set_param(term_block, 'Position', [x_pos, start_y, x_pos+60, start_y+40]);
    
    % Configure as pure resistor (R-only configuration)
    set_param(term_block, 'Resistance', sprintf('%.1f', network_config.physical.Z0));
    set_param(term_block, 'Inductance', '0');
    set_param(term_block, 'Capacitance', 'inf');
    set_param(term_block, 'BranchType', 'R');  % Set to display as R instead of RLC
    
    model_info.blocks.termination = term_block;
    fprintf('  ✓ Added termination load: %.1f Ω\n', network_config.physical.Z0);
    
catch ME
    fprintf('  ✗ Failed to add termination load: %s\n', ME.message);
end

% Position the main ground at the end of the line (like shunt resistor grounds)
try
    ground_block = model_info.blocks.ground;
    x_pos = start_x + network_config.num_segments * block_spacing;
    set_param(ground_block, 'Position', [x_pos, start_y+80, x_pos+50, start_y+110]);
    fprintf('  ✓ Positioned main ground at end of line\n');
catch ME
    fprintf('  ⚠ Could not position main ground: %s\n', ME.message);
end

fprintf('  Network structure complete\n');

end

function connect_network_blocks(model_name, model_info)
%CONNECT_NETWORK_BLOCKS Automatically connect the network blocks

fprintf('  Connecting network blocks...\n');

try
    % Get network configuration
    network_config = model_info.network_config;
    
    % Connect SSTDR system to network
    connect_sstdr_to_network(model_name, model_info);
    
    % Connect transmission line chain
    connect_transmission_line_chain(model_name, network_config, model_info);
    
    % Connect fault elements
    connect_fault_elements(model_name, network_config, model_info);
    
    % Connect termination and ground
    connect_termination_and_ground(model_name, model_info);
    
    fprintf('  ✓ All connections completed\n');
    
catch ME
    fprintf('  ✗ Connection failed: %s\n', ME.message);
    rethrow(ME);
end

end

function connect_sstdr_to_network(model_name, model_info)
%CONNECT_SSTDR_TO_NETWORK Connect SSTDR system to the first transmission line

% Get the first transmission line segment
first_tl = model_info.blocks.network.TL_Segment_1;

if isfield(model_info.blocks, 'sstdr_system')
    % Connect SSTDR subsystem to first transmission line
    sstdr_block = model_info.blocks.sstdr_system;
    sstdr_ground_block = model_info.blocks.sstdr_ground;  % Use dedicated SSTDR ground
    main_ground_block = model_info.blocks.ground;
    
    try
        % Get port handles and debug port structure
        sstdr_ports = get_param(sstdr_block, 'PortHandles');
        tl_ports = get_param(first_tl, 'PortHandles');
        sstdr_ground_ports = get_param(sstdr_ground_block, 'PortHandles');
        main_ground_ports = get_param(main_ground_block, 'PortHandles');
        
        fprintf('    Debug: SSTDR subsystem port structure:\n');
        fprintf('      Outports: %d, LConn: %d, RConn: %d\n', ...
            length(sstdr_ports.Outport), length(sstdr_ports.LConn), length(sstdr_ports.RConn));
        
        % Try different connection strategies based on available ports
        connected = false;
        
        % Strategy 1: If subsystem has electrical connectors (LConn/RConn)
        if length(sstdr_ports.LConn) >= 1 || length(sstdr_ports.RConn) >= 1
            % Handle case where positive port is LConn
            if length(sstdr_ports.LConn) >= 1
                add_line(model_name, sstdr_ports.LConn(1), tl_ports.LConn(1));
                connected = true;
                fprintf('    ✓ Connected SSTDR LConn (positive) to TL\n');
            end
            
            % Handle case where both ports are RConn (positive and negative)
            if length(sstdr_ports.RConn) >= 2
                % First RConn is positive, second is negative
                add_line(model_name, sstdr_ports.RConn(1), tl_ports.LConn(1));
                add_line(model_name, sstdr_ports.RConn(2), sstdr_ground_ports.LConn(1));
                connected = true;
                fprintf('    ✓ Connected SSTDR RConn(1) (positive) to TL\n');
                fprintf('    ✓ Connected SSTDR RConn(2) (negative) to SSTDR ground\n');
            elseif length(sstdr_ports.RConn) >= 1 && length(sstdr_ports.LConn) == 0
                % Only one RConn port (negative) - connect to ground
                add_line(model_name, sstdr_ports.RConn(1), sstdr_ground_ports.LConn(1));
                fprintf('    ✓ Connected SSTDR RConn (negative) to SSTDR ground\n');
            end
            
            % Check for complete connection
            if (length(sstdr_ports.LConn) >= 1 && length(sstdr_ports.RConn) >= 1) || length(sstdr_ports.RConn) >= 2
                fprintf('    ✓ SSTDR subsystem has both positive and negative ports\n');
            end
        end
        
        % Strategy 2: If subsystem has signal outports (for placeholder subsystem)
        if ~connected && length(sstdr_ports.Outport) >= 1
            fprintf('    ⚠ SSTDR subsystem appears to be a placeholder with signal ports\n');
            fprintf('    ⚠ Electrical connection requires proper SSTDR subsystem design\n');
        end
        
        % Always connect TL negative to main ground
        if length(tl_ports.LConn) >= 2
            add_line(model_name, tl_ports.LConn(2), main_ground_ports.LConn(1));
            fprintf('    ✓ Connected TL negative to main ground\n');
        end
        
        if connected
            fprintf('    ✓ SSTDR subsystem successfully connected to network\n');
        else
            fprintf('    ⚠ SSTDR subsystem connection incomplete - check subsystem design\n');
        end
        
    catch ME
        fprintf('    ✗ SSTDR subsystem connection failed: %s\n', ME.message);
        fprintf('    ⚠ Manual connection may be required\n');
    end
    
elseif isfield(model_info.blocks, 'sstdr_source')
    % Connect individual SSTDR blocks
    source_block = model_info.blocks.sstdr_source;
    ground_block = model_info.blocks.ground;
    
    % For Simscape electrical blocks, use port handles
    try
        % Get port handles for proper connection
        source_ports = get_param(source_block, 'PortHandles');
        tl_ports = get_param(first_tl, 'PortHandles');
        ground_ports = get_param(ground_block, 'PortHandles');
        
        % Connect source positive to first TL left positive (port 1)
        add_line(model_name, source_ports.LConn(1), tl_ports.LConn(1));
        
        % Connect source negative to ground
        add_line(model_name, source_ports.RConn(1), ground_ports.LConn(1));
        
        % Connect TL left negative to ground (if second port exists)
        if length(tl_ports.LConn) >= 2
            add_line(model_name, tl_ports.LConn(2), ground_ports.LConn(1));
        end
        
        fprintf('    ✓ Connected voltage source to first transmission line\n');
        
    catch ME
        fprintf('    ✗ Failed to connect source: %s\n', ME.message);
    end
    
    % Connect sensors if they exist
    if isfield(model_info.blocks, 'rx_sensor')
        try
            rx_sensor = model_info.blocks.rx_sensor;
            rx_ports = get_param(rx_sensor, 'PortHandles');
            
            % Connect voltage sensor across the first transmission line
            add_line(model_name, tl_ports.LConn(1), rx_ports.LConn(1));
            if length(tl_ports.LConn) >= 2
                add_line(model_name, tl_ports.LConn(2), rx_ports.RConn(1));
            end
            
            % Connect sensor output to workspace
            if isfield(model_info.blocks, 'rx_output')
                output_ports = get_param(model_info.blocks.rx_output, 'PortHandles');
                add_line(model_name, rx_ports.Outport(1), output_ports.Inport(1));
            end
            
            fprintf('    ✓ Connected voltage sensor\n');
            
        catch ME
            fprintf('    ⚠ Failed to connect voltage sensor: %s\n', ME.message);
        end
    end
    
    fprintf('    ✓ Connected individual SSTDR blocks to network\n');
end

end

function connect_transmission_line_chain(model_name, network_config, model_info)
%CONNECT_TRANSMISSION_LINE_CHAIN Connect transmission line segments in series

for i = 1:(network_config.num_segments - 1)
    current_segment = sprintf('TL_Segment_%d', i);
    next_segment = sprintf('TL_Segment_%d', i + 1);
    
    current_block = model_info.blocks.network.(current_segment);
    next_block = model_info.blocks.network.(next_segment);
    
    try
        % Get port handles for proper connection
        current_ports = get_param(current_block, 'PortHandles');
        next_ports = get_param(next_block, 'PortHandles');
        
        % Check if there's a series fault after this segment
        fault_info = network_config.faults.(sprintf('segment_%d', i));
        
        if strcmp(fault_info.element_type, 'series')
            % Connect through series fault
            fault_name = sprintf('Fault_%d_series', i);
            fault_block = model_info.blocks.network.(fault_name);
            fault_ports = get_param(fault_block, 'PortHandles');
            
            % Current segment right positive -> series fault -> next segment left positive
            add_line(model_name, current_ports.RConn(1), fault_ports.LConn(1));
            add_line(model_name, fault_ports.RConn(1), next_ports.LConn(1));
            
            % Connect negative terminals directly (use second port if available)
            if length(current_ports.RConn) >= 2 && length(next_ports.LConn) >= 2
                add_line(model_name, current_ports.RConn(2), next_ports.LConn(2));
            end
            
            fprintf('    ✓ Connected segments %d-%d through series fault\n', i, i+1);
            
        else
            % Direct connection
            % Current segment right -> next segment left
            add_line(model_name, current_ports.RConn(1), next_ports.LConn(1));
            
            % Connect second wire if available
            if length(current_ports.RConn) >= 2 && length(next_ports.LConn) >= 2
                add_line(model_name, current_ports.RConn(2), next_ports.LConn(2));
            end
            
            fprintf('    ✓ Connected segments %d-%d directly\n', i, i+1);
        end
        
    catch ME
        fprintf('    ✗ Failed to connect segments %d-%d: %s\n', i, i+1, ME.message);
    end
end

end

function connect_fault_elements(model_name, network_config, model_info)
%CONNECT_FAULT_ELEMENTS Connect shunt fault elements to their individual grounds

for i = 1:network_config.num_segments
    fault_info = network_config.faults.(sprintf('segment_%d', i));
    
    if strcmp(fault_info.element_type, 'shunt')
        fault_name = sprintf('Fault_%d_shunt', i);
        segment_name = sprintf('TL_Segment_%d', i);
        ground_name = sprintf('Ground_%d', i);
        
        try
            fault_block = model_info.blocks.network.(fault_name);
            segment_block = model_info.blocks.network.(segment_name);
            shunt_ground = model_info.blocks.network.(ground_name);
            
            fault_ports = get_param(fault_block, 'PortHandles');
            segment_ports = get_param(segment_block, 'PortHandles');
            ground_ports = get_param(shunt_ground, 'PortHandles');
            
            % Connect shunt fault between positive line and its individual ground
            % Connect to the right side (output) of the transmission line segment
            add_line(model_name, segment_ports.RConn(1), fault_ports.LConn(1));
            add_line(model_name, fault_ports.RConn(1), ground_ports.LConn(1));
            
            fprintf('    ✓ Connected shunt fault %d to individual ground\n', i);
            
        catch ME
            fprintf('    ✗ Failed to connect shunt fault %d: %s\n', i, ME.message);
        end
    end
end

end

function connect_termination_and_ground(model_name, model_info)
%CONNECT_TERMINATION_AND_GROUND Connect termination load and ground references

% Get the last transmission line segment
network_config = model_info.network_config;
last_segment_name = sprintf('TL_Segment_%d', network_config.num_segments);
last_segment = model_info.blocks.network.(last_segment_name);

% Connect termination load
if isfield(model_info.blocks, 'termination')
    try
        term_block = model_info.blocks.termination;
        ground_block = model_info.blocks.ground;
        
        % Get port handles
        last_ports = get_param(last_segment, 'PortHandles');
        term_ports = get_param(term_block, 'PortHandles');
        ground_ports = get_param(ground_block, 'PortHandles');
        
        % Connect termination between last segment output and ground
        add_line(model_name, last_ports.RConn(1), term_ports.LConn(1));
        add_line(model_name, term_ports.RConn(1), ground_ports.LConn(1));
        
        % Also connect the negative terminal of last segment to ground
        if length(last_ports.RConn) >= 2
            add_line(model_name, last_ports.RConn(2), ground_ports.LConn(1));
        end
        
        fprintf('    ✓ Connected termination load and ground\n');
        
    catch ME
        fprintf('    ✗ Failed to connect termination: %s\n', ME.message);
    end
end

end

 