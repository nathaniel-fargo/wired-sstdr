function network_config = generate_random_network(varargin)
%GENERATE_RANDOM_NETWORK Generate random network configurations for dataset creation
%
% Usage:
%   config = generate_random_network()  % Default parameters
%   config = generate_random_network('num_segments', 10, 'fault_probability', 0.2)
%
% Optional name-value pairs:
%   'num_segments' - Number of transmission line segments (default: 5-15 random)
%   'fault_probability' - Probability of fault per segment (default: 0.15)
%   'fault_magnitude_range' - [min, max] fault magnitude (default: [0.1, 1.0])
%   'series_bias' - Bias towards series vs shunt faults (default: 0.5, 0=all shunt, 1=all series)
%   'dx' - Segment length (default: 1.0 m)
%   'Z0' - Characteristic impedance (default: 50 Ω)
%   'velocity' - Propagation velocity (default: 2e8 m/s)
%   'seed' - Random seed for reproducibility (default: random)
%
% Output:
%   network_config - Network configuration structure

% Parse optional arguments
p = inputParser;
addParameter(p, 'num_segments', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'fault_probability', 0.15, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'fault_magnitude_range', [0.1, 1.0], @(x) isnumeric(x) && length(x) == 2 && x(1) < x(2));
addParameter(p, 'series_bias', 0.5, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'dx', 1.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'Z0', 50, @(x) isnumeric(x) && x > 0);
addParameter(p, 'velocity', 2e8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'seed', [], @(x) isempty(x) || isnumeric(x));
parse(p, varargin{:});

% Set random seed if provided
if ~isempty(p.Results.seed)
    rng(p.Results.seed);
end

% Determine number of segments
if isempty(p.Results.num_segments)
    num_segments = randi([5, 15]);  % Random between 5 and 15 segments
else
    num_segments = p.Results.num_segments;
end

% Generate load vector
load_vector = zeros(1, num_segments);

for i = 1:num_segments
    if rand() < p.Results.fault_probability
        % Generate a fault
        magnitude = p.Results.fault_magnitude_range(1) + ...
                   rand() * diff(p.Results.fault_magnitude_range);
        
        % Determine if series (positive) or shunt (negative)
        if rand() < p.Results.series_bias
            load_vector(i) = magnitude;  % Series fault
        else
            load_vector(i) = -magnitude; % Shunt fault
        end
    end
    % Otherwise remains 0 (transmission)
end

% Create network configuration
network_config = create_network_config(load_vector, ...
    'dx', p.Results.dx, ...
    'Z0', p.Results.Z0, ...
    'velocity', p.Results.velocity);

% Add generation parameters to metadata
network_config.metadata.generation_params = struct();
network_config.metadata.generation_params.fault_probability = p.Results.fault_probability;
network_config.metadata.generation_params.fault_magnitude_range = p.Results.fault_magnitude_range;
network_config.metadata.generation_params.series_bias = p.Results.series_bias;
if ~isempty(p.Results.seed)
    network_config.metadata.generation_params.seed = p.Results.seed;
end

end

 