function network_config = create_specific_network(network_type, varargin)
%CREATE_SPECIFIC_NETWORK Create specific network patterns for testing
%
% Usage:
%   config = create_specific_network('single_open')
%   config = create_specific_network('multiple_faults', 'num_segments', 10)
%
% Network types:
%   'single_open' - Single open circuit fault in middle
%   'single_short' - Single short circuit fault in middle
%   'multiple_faults' - Multiple random faults
%   'alternating' - Alternating series/shunt faults
%   'gradient' - Increasing fault magnitude
%   'clean' - No faults (pure transmission line)

% Parse arguments
p = inputParser;
addRequired(p, 'network_type', @ischar);
addParameter(p, 'num_segments', 8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fault_magnitude', 0.8, @(x) isnumeric(x) && x > 0 && x <= 1);
parse(p, network_type, varargin{:});

num_segments = p.Results.num_segments;
fault_mag = p.Results.fault_magnitude;

switch lower(network_type)
    case 'single_open'
        load_vector = zeros(1, num_segments);
        mid_point = ceil(num_segments / 2);
        load_vector(mid_point) = fault_mag;
        
    case 'single_short'
        load_vector = zeros(1, num_segments);
        mid_point = ceil(num_segments / 2);
        load_vector(mid_point) = -fault_mag;
        
    case 'multiple_faults'
        load_vector = zeros(1, num_segments);
        fault_positions = randperm(num_segments, min(3, num_segments));
        for pos = fault_positions
            if rand() > 0.5
                load_vector(pos) = fault_mag * (0.5 + 0.5 * rand());
            else
                load_vector(pos) = -fault_mag * (0.5 + 0.5 * rand());
            end
        end
        
    case 'alternating'
        load_vector = zeros(1, num_segments);
        for i = 2:2:num_segments
            if mod(i/2, 2) == 1
                load_vector(i) = fault_mag;  % Series
            else
                load_vector(i) = -fault_mag; % Shunt
            end
        end
        
    case 'gradient'
        load_vector = zeros(1, num_segments);
        for i = 1:num_segments
            if i > num_segments/3 && i < 2*num_segments/3
                magnitude = fault_mag * (i - num_segments/3) / (num_segments/3);
                load_vector(i) = magnitude;
            end
        end
        
    case 'clean'
        load_vector = zeros(1, num_segments);
        
    otherwise
        error('Unknown network type: %s', network_type);
end

% Create configuration (filter out our local parameters)
config_args = {};
for i = 1:2:length(varargin)
    if ~strcmp(varargin{i}, 'num_segments') && ~strcmp(varargin{i}, 'fault_magnitude')
        config_args{end+1} = varargin{i};
        config_args{end+1} = varargin{i+1};
    end
end

network_config = create_network_config(load_vector, config_args{:});
network_config.metadata.network_type = network_type;

end 