function [R_value, element_type] = map_load_to_resistance(load_value, varargin)
%MAP_LOAD_TO_RESISTANCE Convert load value to resistance for Simscape elements
%
% Usage:
%   [R, type] = map_load_to_resistance(0.5)     % Series resistance
%   [R, type] = map_load_to_resistance(-0.3)    % Shunt resistance
%   [R, type] = map_load_to_resistance(0, 'Z0', 75)  % Custom characteristic impedance
%
% Input:
%   load_value - Scalar in range [-1, 1]
%                0: transmission (no extra element)
%                >0: series fault (higher = more open)
%                <0: shunt fault (more negative = more short)
%   varargin   - Optional name-value pairs:
%                'Z0' - Characteristic impedance for scaling (default: 50 Ω)
%                'series_range' - [R_min, R_max] for series faults (default: [1e3, 1e6])
%                'shunt_range' - [R_min, R_max] for shunt faults (default: [1, 1e3])
%                'mapping' - 'linear' or 'log' (default: 'log')
%
% Output:
%   R_value      - Resistance value in Ohms
%   element_type - 'none', 'series', or 'shunt'

% Parse optional arguments
p = inputParser;
addParameter(p, 'Z0', 50, @(x) isnumeric(x) && x > 0);
addParameter(p, 'series_range', [1e3, 1e6], @(x) isnumeric(x) && length(x) == 2 && x(1) < x(2));
addParameter(p, 'shunt_range', [1, 1e3], @(x) isnumeric(x) && length(x) == 2 && x(1) < x(2));
addParameter(p, 'mapping', 'log', @(x) ismember(x, {'linear', 'log'}));
parse(p, varargin{:});

% Validate input
if load_value < -1 || load_value > 1
    error('map_load_to_resistance:OutOfRange', ...
        'Load value must be in range [-1, 1], got %.3f', load_value);
end

% Determine element type and resistance value
if abs(load_value) < 1e-6
    % Transmission line segment (no extra element)
    R_value = NaN;
    element_type = 'none';
    
elseif load_value > 0
    % Series fault (positive load = more open)
    element_type = 'series';
    R_min = p.Results.series_range(1);
    R_max = p.Results.series_range(2);
    
    if strcmp(p.Results.mapping, 'linear')
        R_value = R_min + load_value * (R_max - R_min);
    else % logarithmic
        log_min = log10(R_min);
        log_max = log10(R_max);
        R_value = 10^(log_min + load_value * (log_max - log_min));
    end
    
else
    % Shunt fault (negative load = more short)
    element_type = 'shunt';
    R_min = p.Results.shunt_range(1);
    R_max = p.Results.shunt_range(2);
    
    % Use absolute value and invert mapping (more negative = lower resistance)
    abs_load = abs(load_value);
    
    if strcmp(p.Results.mapping, 'linear')
        R_value = R_max - abs_load * (R_max - R_min);
    else % logarithmic
        log_min = log10(R_min);
        log_max = log10(R_max);
        R_value = 10^(log_max - abs_load * (log_max - log_min));
    end
end

% Ensure minimum resistance for numerical stability
if ~isnan(R_value) && R_value < 0.1
    R_value = 0.1;
end

end 