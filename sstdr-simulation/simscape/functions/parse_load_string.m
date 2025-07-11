function [load_vector, metadata] = parse_load_string(load_string, varargin)
%PARSE_LOAD_STRING Convert load string notation to numeric vector
%
% Usage:
%   [loads, meta] = parse_load_string('(0)(0.3)(-0.5)(1)')
%   [loads, meta] = parse_load_string('(0)(0.3)(-0.5)(1)', 'dx', 1.0, 'Z0', 50)
%
% Input:
%   load_string - String in format "(val1)(val2)...(valN)" where each val ∈ [-1,1]
%   varargin    - Optional name-value pairs:
%                 'dx' - Physical distance per segment (default: 1.0 m)
%                 'Z0' - Characteristic impedance (default: 50 Ω)
%                 'velocity' - Propagation velocity (default: 2e8 m/s)
%
% Output:
%   load_vector - Numeric vector of load values
%   metadata    - Struct with physical parameters and validation info

% Parse optional arguments
p = inputParser;
addParameter(p, 'dx', 1.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'Z0', 50, @(x) isnumeric(x) && x > 0);
addParameter(p, 'velocity', 2e8, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});

% Extract load values using regex
pattern = '\(([+-]?\d*\.?\d+)\)';
matches = regexp(load_string, pattern, 'tokens');

if isempty(matches)
    error('parse_load_string:InvalidFormat', ...
        'No valid load values found. Expected format: (val1)(val2)...(valN)');
end

% Convert to numeric vector
load_vector = zeros(1, length(matches));
for i = 1:length(matches)
    val = str2double(matches{i}{1});
    if isnan(val)
        error('parse_load_string:InvalidValue', ...
            'Invalid numeric value at position %d: %s', i, matches{i}{1});
    end
    if val < -1 || val > 1
        error('parse_load_string:OutOfRange', ...
            'Load value at position %d is out of range [-1,1]: %.3f', i, val);
    end
    load_vector(i) = val;
end

% Create metadata structure
metadata = struct();
metadata.num_segments = length(load_vector);
metadata.dx = p.Results.dx;
metadata.Z0 = p.Results.Z0;
metadata.velocity = p.Results.velocity;
metadata.total_length = metadata.num_segments * metadata.dx;
metadata.one_way_delay = metadata.total_length / metadata.velocity;
metadata.original_string = load_string;

% Validation summary
metadata.validation = struct();
metadata.validation.num_opens = sum(load_vector > 0.5);
metadata.validation.num_shorts = sum(load_vector < -0.5);
metadata.validation.num_transmission = sum(abs(load_vector) < 0.1);
metadata.validation.max_load = max(abs(load_vector));

fprintf('Parsed load string: %d segments, %.1f m total length\n', ...
    metadata.num_segments, metadata.total_length);
fprintf('  Opens: %d, Shorts: %d, Transmission: %d\n', ...
    metadata.validation.num_opens, metadata.validation.num_shorts, ...
    metadata.validation.num_transmission);

end 