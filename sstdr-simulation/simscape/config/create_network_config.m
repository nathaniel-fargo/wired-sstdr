function network_config = create_network_config(load_vector, varargin)
%CREATE_NETWORK_CONFIG Create a network configuration structure
%
% Usage:
%   config = create_network_config([0, 0.5, -0.3, 0])
%   config = create_network_config([0, 0.5, -0.3], 'dx', 0.5, 'Z0', 75)
%
% Input:
%   load_vector - Numeric vector of load values ∈ [-1, 1]
%                 0: transmission (no fault)
%                 >0: series fault (higher = more open)
%                 <0: shunt fault (more negative = more short)
%   varargin    - Optional name-value pairs:
%                 'dx' - Physical distance per segment (default: 1.0 m)
%                 'Z0' - Characteristic impedance (default: 50 Ω)
%                 'velocity' - Propagation velocity (default: 2e8 m/s)
%                 'termination_impedance' - Termination load impedance (default: Z0)
%                 'name' - Network name (default: auto-generated)
%
% Output:
%   network_config - Structure containing all network parameters

% Parse optional arguments
p = inputParser;
addRequired(p, 'load_vector', @(x) isnumeric(x) && isvector(x));
addParameter(p, 'dx', 1.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'Z0', 50, @(x) isnumeric(x) && x > 0);
addParameter(p, 'velocity', 2e8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'termination_impedance', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'name', '', @ischar);
parse(p, load_vector, varargin{:});

% Validate load vector
if any(load_vector < -1 | load_vector > 1)
    error('create_network_config:InvalidRange', ...
        'All load values must be in range [-1, 1]');
end

% Ensure row vector
load_vector = load_vector(:)';

% Generate name if not provided
if isempty(p.Results.name)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    network_name = sprintf('network_%s', timestamp);
else
    network_name = p.Results.name;
end

% Create network configuration structure
network_config = struct();

% Basic network parameters
network_config.name = network_name;
network_config.type = 'line';  % For future extension to branching
network_config.num_segments = length(load_vector);
network_config.load_vector = load_vector;

% Physical parameters
network_config.physical = struct();
network_config.physical.dx = p.Results.dx;
network_config.physical.Z0 = p.Results.Z0;
network_config.physical.velocity = p.Results.velocity;
network_config.physical.total_length = network_config.num_segments * p.Results.dx;
network_config.physical.one_way_delay = network_config.physical.total_length / p.Results.velocity;

% Handle termination impedance from final load vector value
final_load_value = load_vector(end);
[term_resistance, ~] = map_load_to_resistance(final_load_value, 'Z0', p.Results.Z0);

% Set termination impedance based on final load vector value
if isempty(p.Results.termination_impedance)
    % Use final load vector value to determine termination impedance
    if abs(final_load_value) < 1e-6
        % Load value ~0: matched termination (Z0)
        network_config.physical.termination_impedance = p.Results.Z0;
    elseif final_load_value > 0
        % Positive load value: high impedance (open-like termination)
        network_config.physical.termination_impedance = term_resistance;
    else
        % Negative load value: low impedance (short-like termination)
        network_config.physical.termination_impedance = term_resistance;
    end
else
    % Use explicitly specified termination impedance
    network_config.physical.termination_impedance = p.Results.termination_impedance;
end

% Analysis of network characteristics (exclude termination from fault count)
load_vector_faults = load_vector(1:end-1);  % Exclude final termination value
network_config.analysis = struct();
network_config.analysis.num_faults = sum(abs(load_vector_faults) > 1e-6);
network_config.analysis.num_opens = sum(load_vector_faults > 0.1);
network_config.analysis.num_shorts = sum(load_vector_faults < -0.1);
network_config.analysis.num_transmission = sum(abs(load_vector_faults) < 0.1);
network_config.analysis.max_fault_magnitude = max(abs(load_vector_faults));
network_config.analysis.fault_positions = find(abs(load_vector_faults) > 1e-6);
network_config.analysis.termination_load_value = final_load_value;
if final_load_value > 0.1
    network_config.analysis.termination_type = 'open';
elseif final_load_value < -0.1
    network_config.analysis.termination_type = 'short';
else
    network_config.analysis.termination_type = 'matched';
end

% Fault details for each segment (excluding the final termination segment)
network_config.faults = struct();
for i = 1:(network_config.num_segments - 1)  % Exclude final segment
    segment_name = sprintf('segment_%d', i);
    
    [R_value, element_type] = map_load_to_resistance(load_vector(i));
    
    network_config.faults.(segment_name) = struct();
    network_config.faults.(segment_name).load_value = load_vector(i);
    network_config.faults.(segment_name).element_type = element_type;
    network_config.faults.(segment_name).resistance = R_value;
    network_config.faults.(segment_name).position_m = i * p.Results.dx;
    network_config.faults.(segment_name).delay_s = network_config.faults.(segment_name).position_m / p.Results.velocity;
end

% Store termination information separately
network_config.termination = struct();
network_config.termination.load_value = final_load_value;
network_config.termination.resistance = network_config.physical.termination_impedance;
network_config.termination.position_m = network_config.num_segments * p.Results.dx;
network_config.termination.delay_s = network_config.termination.position_m / p.Results.velocity;

% Metadata
network_config.metadata = struct();
network_config.metadata.created = datestr(now);
network_config.metadata.version = '1.0';
network_config.metadata.description = sprintf('%d-segment 1D transmission line network', network_config.num_segments);

% Summary display
fprintf('Created network configuration: %s\n', network_config.name);
fprintf('  Segments: %d, Length: %.1f m\n', network_config.num_segments, network_config.physical.total_length);
fprintf('  Faults: %d (%d opens, %d shorts)\n', ...
    network_config.analysis.num_faults, network_config.analysis.num_opens, network_config.analysis.num_shorts);
fprintf('  Max fault magnitude: %.2f\n', network_config.analysis.max_fault_magnitude);

end 